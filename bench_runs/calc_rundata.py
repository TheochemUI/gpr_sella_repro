import click
from pathlib import Path
import pandas as pd
import ase.io
import dataclasses


from chemparseplot.parse.eon.saddle_search import parse_eon_saddle
from chemparseplot.parse.eon.minimization import min_e_result
from chemparseplot.parse.sella.saddle_search import (
    parse_sella_saddle,
    _no_ghost,
    _get_ghosts,
)
from chemparseplot.basetypes import SpinID
from chemparseplot.analyze.use_ira import calculate_rmsd

SINGLET_IDS = [f"{x:003d}" for x in range(265)]
DOUBLET_IDS = [f"{x:003d}" for x in range(235)]
SPINS = ["singlets", "doublets"]


def add_min(resdat, mindir, init, fin):
    mindir = Path(mindir)
    if (mindir / "min.con").exists():
        mincon = ase.io.read(mindir / "min.con")
        resdat |= [
            ("min_energy", min_e_result(mindir)),
            ("rmsd_min_sad", calculate_rmsd(mincon, fin, 10)),
            ("rmsd_min_init", calculate_rmsd(mincon, init, 10)),
        ]
    return resdat


def get_eon_result(basepath: Path, spindat: SpinID, mindir: Path):
    rpath = basepath / "results.dat"
    result_data = dataclasses.asdict(parse_eon_saddle(basepath, spindat))
    if rpath.exists():
        init = ase.io.read(basepath / "pos.con")
        fin = ase.io.read(basepath / "saddle.con")
        result_data |= [
            ("rmsd_init_saddle", calculate_rmsd(init, fin, 10)),
            ("nghost", 0),
        ]
        if mindir:
            return add_min(result_data, mindir, init, fin)
    else:
        click.echo(f"No result in {basepath}", err=True)
    return result_data


def get_sella_result(basepath: Path, spindat: SpinID, mindir: Path):
    rtraj = list(basepath.glob("*.traj"))
    if not len(rtraj) == 1:
        click.echo(f"No result in {basepath}", err=True)
    result_data = dataclasses.asdict(parse_sella_saddle(basepath, spindat))
    ctraj = rtraj[0]
    if ctraj.exists():
        init = _no_ghost(ase.io.read(ctraj, 0))
        fin = _no_ghost(ase.io.read(ctraj, -1))
        result_data |= [
            ("rmsd_init_saddle", calculate_rmsd(init, fin, 10)),
            ("nghost", _get_ghosts(ctraj)),
        ]
        if mindir:
            return add_min(result_data, mindir, init, fin)
    return result_data


@click.command()
@click.option(
    "--base_path",
    default=".",
    help="Base path for runs. Usually the directory containing the 'runs' folder.",
)
@click.option(
    "--code",
    default="eon",
    type=click.Choice(["eon", "sella"]),
    help="Code used for generating runs.",
)
@click.option("--scf_thresh", default="1e8m", help="SCF convergence threshold.")
@click.option(
    "--run_dir",
    default="default_run",
    help="Name of the run directory, with the method, like idimer/blah.",
)
@click.option(
    "--min_dir",
    default=None,
    help="Name of the directory which has minimization data, like min/hpc_min.",
)
@click.option("--output", default="data_impd.csv", help="Output CSV file name.")
def cli(base_path, scf_thresh, run_dir, min_dir, output, code):
    """Parses saddle point search results and saves them to a CSV file."""

    base_path = Path(base_path)
    all_results = []

    for spin in SPINS:
        mol_ids = SINGLET_IDS if spin == "singlets" else DOUBLET_IDS
        for mol_id in mol_ids:
            click.echo(f"Processing: {scf_thresh}/{run_dir}/{spin}/{mol_id}")
            results_path = base_path / f"runs/{scf_thresh}/{run_dir}/{spin}/{mol_id}"
            if min_dir:
                min_path = f"{Path(min_dir).absolute()}/{spin}/{mol_id}"
            if code.lower() == "eon":
                all_results.append(
                    get_eon_result(results_path, SpinID(mol_id, spin), min_path)
                )
            elif code.lower() == "sella":
                all_results.append(
                    get_sella_result(results_path, SpinID(mol_id, spin), min_path)
                )
            else:
                click.echo(
                    f"Skipping as no result in {results_path} for {code}", err=True
                )

    # Generate dataframe
    df = pd.DataFrame(all_results)
    if not df.empty:
        df = df.sort_values(["mol_id", "spin"], ignore_index=True)
    df.to_csv(output, index=False, encoding="utf-8")
    click.echo(f"Data saved to {output}")


if __name__ == "__main__":
    cli()
