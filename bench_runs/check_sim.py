from pathlib import Path


import ase.io
import click
import numpy as np
import pandas as pd
import ira_mod
from chemparseplot.analyze.dist import analyze_structure

SINGLET_IDS = [f"{x:003d}" for x in range(265)]
DOUBLET_IDS = [f"{x:003d}" for x in range(235)]
SPINS = ["singlets", "doublets"]


def is_sim_pair(atm1, atm2, e_tol=0.01, hd_tol=1, k_factor=1.8):
    ira = ira_mod.IRA()
    # XXX: the atoms need the potential energy attached, i.e these are from a trajectory
    if len(atm1)!=len(atm2):
        return False
    if np.all(atm1.symbols == atm2.symbols):
        if (
            abs(abs(atm1.get_potential_energy()) - abs(atm2.get_potential_energy()))
            < e_tol
        ):
            hd = ira.match(
                len(atm1),
                atm1.get_atomic_numbers(),
                atm1.get_positions(),
                len(atm2),
                atm2.get_atomic_numbers(),
                atm2.get_positions(),
                k_factor,
            )[-1]
            if hd < hd_tol:
                return True
            else:
                return False
        else:
            return False
    else:
        return False

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
            if code.lower() == "sella":
              traj_files = list(results_path.glob("*.traj"))
            else:
              traj_files = list(results_path.glob(f"*_{code.lower()}.traj")) #For this use case.
            if not traj_files:
                click.echo(f"No .traj files found in {results_path}", err=True)
                continue
            # Select the *first* .traj file, and error if it doesn't.
            traj_file = traj_files[0]

            try:
                # Read *only* the first frame using index=0, else -1 for saddle
                atoms = ase.io.read(str(traj_file), index=0)
                atoms_dict[(mol_id, spin)] = atoms
                _,_,ifrags,hdists_i,cdists_i = analyze_structure(atoms)
                # _,_,sfrags,hdists_s = analyze_structure(ase.io.read(str(traj_file), index=-1))
                print(f"{spin}, {mol_id} :: {len(ifrags)} {np.max(hdists_i):02f} {cdists_i}")
            except Exception as e:
                click.echo(f"Error reading trajectory {traj_file}: {e}", err=True)
                # Don't return, continue to the next.
                continue
    return atoms_dict


@click.command()
@click.option("--base_path", default=".", type=click.Path(exists=True, file_okay=False, dir_okay=True, path_type=Path),
              help="Base path for runs (directory containing 'runs').")
@click.option("--scf_thresh", default="1e8m", help="SCF convergence threshold.")
@click.option("--run_dir", default="default_run", help="Name of the run directory.")
@click.option("--output", default="duplicate_systems.csv", help="Output CSV file name.")
@click.option("--tolerance", "-tol", default=1e-5, type=float,
              help="Tolerance for comparing the energy of atoms.")
@click.option("--code", default="sella", type=str, help="Code for parsing the file")
def find_duplicates(base_path: Path, scf_thresh: str, run_dir: str, output: str, tolerance: float, code:str):
    """Finds and reports duplicate structures."""

    click.echo(f"Loading trajectories from: {base_path / 'runs' / scf_thresh / run_dir}")
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
        for j in range(i + 1, len(atoms_list)):  # Start j from i+1 to avoid redundant comparisons
            (mol_id1, spin1), atoms1 = atoms_list[i]
            (mol_id2, spin2), atoms2 = atoms_list[j]

            pair_id = tuple(sorted(((mol_id1, spin1), (mol_id2, spin2)))) #make the id
            if pair_id in processed_pairs:
                continue #skip
            processed_pairs.add(pair_id)


            if is_sim_pair(atoms1, atoms2, tolerance):
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
        click.echo(f"Found {len(duplicates)} duplicate systems. Results saved to {output}")
    else:
        click.echo("No duplicate systems found.")


if __name__ == "__main__":
    find_duplicates()
