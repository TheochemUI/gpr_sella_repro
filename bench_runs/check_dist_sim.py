# eg. python check_dist_sim.py --base_path eon --run_dir "min/hpc_min" --tolerance 0.2 --code eon
from pathlib import Path


import ase.io
import click
import numpy as np
import pandas as pd
from chemparseplot.analyze.dist import analyze_structure
from chemparseplot.analyze.use_ira import is_ira_pair

SINGLET_IDS = [f"{x:003d}" for x in range(265)]
DOUBLET_IDS = [f"{x:003d}" for x in range(235)]
SPINS = ["singlets", "doublets"]


def load_trajectories(base_path: Path, scf_thresh: str, run_dir: str, code: str):
    """Loads trajectories and returns a dictionary of Atoms objects.

    Args:
        base_path: Base path for the runs.
        scf_thresh: SCF threshold directory name.
        run_dir: Run directory name.

    Returns:
        A dictionary where keys are (mol_id, spin) tuples and values are
        the corresponding ASE Atoms objects (last frame of the trajectory).
        Returns an empty dictionary if any error occurs.
    """
    atoms_dict = {}
    for spin in SPINS:
        mol_ids = SINGLET_IDS if spin == "singlets" else DOUBLET_IDS
        for mol_id in mol_ids:
            results_path = base_path / "runs" / scf_thresh / run_dir / spin / mol_id
            if not results_path.exists():
                click.echo(f"Skipping missing directory: {results_path}", err=True)
                continue
            if code.lower() == "eon":
                mincon_files = list(results_path.glob("min.con"))
            if not mincon_files:
                click.echo(f"No .con files found in {results_path}", err=True)
                continue
            minf = mincon_files[0]

            try:
                # Read *only* the first frame using index=0, else -1 for saddle
                atoms = ase.io.read(minf)
                atoms_dict[(mol_id, spin)] = atoms
                _, _, ifrags, hdists_i, cdists_i = analyze_structure(atoms)
                print(
                    f"{spin}, {mol_id} :: {len(ifrags)} {np.max(hdists_i):02f} {cdists_i}"
                )
            except Exception as e:
                click.echo(f"Error reading trajectory {minf}: {e}", err=True)
                # Don't return, continue to the next.
                continue
    return atoms_dict


@click.command()
@click.option(
    "--base_path",
    default=".",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
    help="Base path for runs (directory containing 'runs').",
)
@click.option("--scf_thresh", default="1e8m", help="SCF convergence threshold.")
@click.option("--run_dir", default="default_run", help="Name of the run directory.")
@click.option("--output", default="duplicate_systems.csv", help="Output CSV file name.")
@click.option(
    "--tolerance",
    "-tol",
    default=1e-5,
    type=float,
    help="Tolerance for comparing the energy of atoms.",
)
@click.option("--code", default="sella", type=str, help="Code for parsing the file")
def find_duplicates(
    base_path: Path,
    scf_thresh: str,
    run_dir: str,
    output: str,
    tolerance: float,
    code: str,
):
    """Finds and reports duplicate structures."""

    click.echo(
        f"Loading trajectories from: {base_path / 'runs' / scf_thresh / run_dir}"
    )
    atoms_dict = load_trajectories(base_path, scf_thresh, run_dir, code)

    if not atoms_dict:
        click.echo("No valid trajectories loaded. Exiting.", err=True)
        return

    click.echo("Comparing structures...")
    duplicates = []
    processed_pairs = set()  # Keep track of pairs we've already compared

    # Convert dictionary to list for easier iteration
    atoms_list = list(atoms_dict.items())

    for i in range(len(atoms_list)):
        for j in range(
            i + 1, len(atoms_list)
        ):  # Start j from i+1 to avoid redundant comparisons
            (mol_id1, spin1), atoms1 = atoms_list[i]
            (mol_id2, spin2), atoms2 = atoms_list[j]

            pair_id = tuple(sorted(((mol_id1, spin1), (mol_id2, spin2))))  # make the id
            if pair_id in processed_pairs:
                continue  # skip
            processed_pairs.add(pair_id)

            if is_ira_pair(atoms1, atoms2, tolerance):
                duplicates.append(
                    {
                        "mol_id1": mol_id1,
                        "spin1": spin1,
                        "mol_id2": mol_id2,
                        "spin2": spin2,
                    }
                )

    if duplicates:
        df = pd.DataFrame(duplicates)
        df.to_csv(output, index=False)
        click.echo(
            f"Found {len(duplicates)} duplicate systems. Results saved to {output}"
        )
    else:
        click.echo("No duplicate systems found.")


if __name__ == "__main__":
    find_duplicates()
