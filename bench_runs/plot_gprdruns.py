#!/usr/bin/env python3
from pathlib import Path

from chemparseplot.parse.eon.gprd import (
    create_geom_traj_from_hdf5,
    create_full_traj_from_hdf5,
    create_nwchem_trajectory,
)
import ase.io as aseio

import click


@click.group()
def cli():
    """A CLI for creating ASE trajectories from HDF5 files."""
    pass


# TODO: Maybe this ought to live within the chemparseplot code


@cli.command()
@click.argument(
    "hdf5_file", type=click.Path(exists=True, dir_okay=False, path_type=Path)
)
@click.option(
    "--outer_loop_group", default="outer_loop", help="Name of the HDF5 group."
)
@click.option(
    "--initial_structure",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to initial structure. Defaults to 'pos.con' in HDF5 dir.",
)
@click.option(
    "--output_dir",
    "-odir",
    type=click.Path(file_okay=False, path_type=Path),
    help="Output directory. Defaults to the HDF5 file directory.",
)
@click.option(
    "--output_name",
    "-oname",
    type=str,
    help="Base name for output trajectory. Defaults to HDF5 file stem.",
)
@click.option(
    "--traj_type",
    "-tt",
    type=click.Choice(["geom", "full"]),
    help="Base name for output trajectory. Defaults to HDF5 file stem.",
)
def hdf5_to_traj(
    hdf5_file, outer_loop_group, initial_structure, output_dir, output_name, traj_type
):
    """Creates an ASE trajectory file from an HDF5 file."""

    # --- 1. Handle Default Paths ---
    if initial_structure is None:
        initial_structure = hdf5_file.parent / "pos.con"
    if not initial_structure.exists():
        raise click.ClickException(
            f"Initial structure file not found: {initial_structure}"
        )

    if output_dir is None:
        output_dir = hdf5_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    if output_name is None:
        output_name = hdf5_file.stem
    output_traj_file = output_dir / (f"{ output_name.replace('.traj', '') }.traj")

    # --- 2. Create Trajectory ---
    click.echo(f"Creating trajectory from HDF5: {hdf5_file} -> {output_traj_file}")
    if traj_type == "geom":
        create_geom_traj_from_hdf5(
            hdf5_file, output_traj_file, str(initial_structure), outer_loop_group
        )
    elif traj_type == "full":
        create_full_traj_from_hdf5(
            hdf5_file, output_traj_file, str(initial_structure), outer_loop_group
        )

    click.echo("Trajectory creation complete.")


@cli.command()
@click.option(
    "--hdf5_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the HDF5 file.",
)
@click.option(
    "--output_traj",
    required=True,
    type=click.Path(dir_okay=False, path_type=Path),
    help="Path to the output trajectory file.",
)
@click.option(
    "--output_dir",
    "-odir",
    type=click.Path(file_okay=False, path_type=Path),
    help="Output directory. Defaults to the HDF5 file directory.",
)
@click.option(
    "--initial_structure",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the initial structure file.",
)
@click.option(
    "--outer_loop_group",
    default="outer_loop",
    help="Name of the HDF5 group containing outer loop data. Defaults to 'outer_loop'.",
)
@click.option(
    "--multiplicity",
    default=1,
    type=click.Choice([1, 2]),
    help="Multiplicity for NWChem calculation",
)
def create_nwchem(
    hdf5_file,
    output_traj,
    output_dir,
    initial_structure,
    outer_loop_group,
    multiplicity,
):
    """Creates an ASE trajectory file and calculates energy/forces with NWChem."""
    # --- 1. Handle Default Paths ---
    if initial_structure is None:
        initial_structure = hdf5_file.parent / "pos.con"
    if not initial_structure.exists():
        raise click.ClickException(
            f"Initial structure file not found: {initial_structure}"
        )

    if output_dir is None:
        output_dir = hdf5_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    if output_traj is None:
        output_traj = hdf5_file.stem
    output_traj_file = output_dir / (f"{ output_traj.replace('.traj', '') }.traj")

    try:
        initial_atoms = aseio.read(initial_structure)
    except FileNotFoundError:
        print(f"Error: Could not read initial structure from '{initial_structure}'")
        return
    except Exception as e:
        print(f"An error occurred while reading initial structure: {e}")
        return

    create_nwchem_trajectory(
        initial_atoms,
        str(hdf5_file),
        str(output_traj_file),
        multiplicity,
        outer_loop_group,
    )


if __name__ == "__main__":
    cli()
