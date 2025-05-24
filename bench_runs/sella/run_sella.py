#!/usr/bin/env python3

import os
from io import StringIO
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import shutil
import logging

import ase.io
from ase.io import read as aseread
from ase.io.trajectory import Trajectory as ASETraj
from ase.calculators.nwchem import NWChem
from ase.calculators.socketio import SocketIOCalculator
from sella import Sella

from rgpycrumbs._aux import getstrform, switchdir

import eon.akmc
from eon.explorer import MinModeExplorer
from eon.config import ConfigClass as eonConf
import eon.fileio as eio

import click


def setup_logger():
    """Sets up the logger."""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    return logger


def prepare_min_con(run_dir, spin, mol_idx_str, min_dir):
    """Prepares the pos.con from a min.con file."""
    atm = aseread(Path(f"{min_dir}/{spin}/{mol_idx_str}/min.con"))
    atm.set_cell([25, 25, 25])
    atm.center()
    ase.io.write(getstrform(Path(run_dir) / "pos.con"), atm)


def setup_and_run_initial_displacement(run_dir, logger):
    """Sets up and runs the initial displacement."""
    with switchdir(run_dir):
        econf = eonConf()
        econf.init("base_config.ini")

        kT = econf.main_temperature / 11604.5  # in eV
        states = eon.akmc.get_statelist(kT, econf)
        start_state_num = 0
        current_state = states.get_state(start_state_num)
        previous_state = current_state
        explore_state = current_state
        state_explorer = MinModeExplorer(
            states, previous_state, explore_state, superbasin=None, config=econf
        )
        displacement, mode, disp_type = state_explorer.generate_displacement()
        eio.savecon(Path.cwd() / "displacement.con", displacement)
        eio.save_mode(Path.cwd() / "direction.dat", mode)

        # Delete eOn's generated files that we don't need here
        shutil.rmtree(Path.cwd() / "jobs")
        shutil.rmtree(Path.cwd() / "states")

        logger.info(f"Initial displacement setup complete in {run_dir}")


def prepare_atom(spin, mol_idx_str, datzip):
    """Prepares the pos.con file."""
    with ZipFile(datzip, "r") as zdat:
        with zdat.open(f"sella_si/{spin}/{mol_idx_str}.xyz", "r") as atmdat:
            atoms = aseread(StringIO(atmdat.read().decode()), format="xyz")
            return atoms


def convert_scfthresh(scfthresh_input):
    """Converts scfthresh input (e.g., 1e8m) to config value (e.g., 1e-8).

    Args:
        scfthresh_input: The input string (e.g., "1e8m", "1e-5").

    Returns:
        The converted scfthresh value (e.g., "1e-8", "1e-5").

    Raises:
        ValueError: If the input format is invalid.
    """
    if scfthresh_input.endswith("m"):
        try:
            base, exponent = scfthresh_input[:-1].split("e")
            return f"{base}e{int(exponent) * -1}"
        except ValueError:
            raise ValueError(f"Invalid scfthresh input format: {scfthresh_input}")
    elif "e" in scfthresh_input and not scfthresh_input.endswith("m"):
        return scfthresh_input
    else:
        raise ValueError(f"Invalid scfthresh input format: {scfthresh_input}")


@click.command()
@click.argument("mol_idx", type=int)
@click.option(
    "--spin",
    type=click.Choice(["singlets", "doublets"]),
    default="singlets",
    help="Multiplicity of the molecule (singlets or doublets)",
)
@click.option(
    "--scfthresh",
    type=click.Choice(["1e8m", "1e5m"]),
    default="1e8m",
    help="SCF Threshold",
)
@click.option(
    "--rdir",
    type=str,
    default="smol_hpc",
    help="Run directory name within runs/{scfthresh}/",
)
@click.option(
    "--datzip",
    type=str,
    default="/tmp/nwchem_scratch",
    help="Zip file with sella_si/{spin}/{mol_idx_str}.xyz",
)
@click.option(
    "--min_dir",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, path_type=Path
    ),
    default=None,
    help="Directory containing a min.con file to use instead of generating pos.con.",
)
def main(mol_idx, spin, scfthresh, rdir, datzip, min_dir):
    """GPRD benchmarks on the Sella test set.

    MOL_IDX is the index of the molecule to optimize (e.g., 011).
    """
    logger = setup_logger()

    mult = 1 if spin == "singlets" else 2
    mol_idx_str = f"{mol_idx:03d}"

    # Full path to your nwchem executable
    nwchem_path = os.environ["NWCHEM_COMMAND"]

    # Memory for NWChem to allocate. You probably don't need to change this.
    memory = "2 gb"

    label = f"{spin}_{mol_idx_str}"

    nwchem_kwargs = dict(
        label=label,
        directory=".",
        set={"geom:dont_verify": True},
        command=(
            f"mpirun -n {os.environ['NWCHEM_MPI']}"
            f" {nwchem_path} PREFIX.nwi > PREFIX.nwo"
        ),
        memory=memory,
        scf=dict(
            nopen=mult - 1,
            # 1e-5 takes longer! (almost twice as much)
            thresh=float(convert_scfthresh(scfthresh)),
            maxiter=200,
        ),
        basis="3-21G",
        task="optimize",
        driver=dict(
            socket=f"unix {label}",
        ),
    )

    if mult == 2:
        nwchem_kwargs["scf"]["uhf"] = None

    nwchem = NWChem(**nwchem_kwargs)

    logger.info(f"Running mult({mult}) :: {mol_idx_str} in {Path.cwd()}")
    if min_dir:
        prepare_min_con(
            ".",
            spin,
            mol_idx_str,
            min_dir,
        )
        setup_and_run_initial_displacement(".", logger)
        atm = aseread(Path("displacement.con"))
        # Now we construct the same initial starting point as the dimer methods
        eigenvector_init = np.loadtxt("direction.dat")
        atm.positions = atm.positions + 0.01 * eigenvector_init
    else:
        atm = prepare_atom(spin, mol_idx_str, datzip)

    with SocketIOCalculator(nwchem, unixsocket=label) as calc:
        atm.calc = calc
        opt = Sella(
            atm,
            internal=True,
            logfile=f"{label}.log",
            trajectory=f"{label}.traj",
        )
        opt.run(fmax=0.01, steps=200)
        npes = len(ASETraj(f"{label}.traj"))
        Path("npes.txt").write_text(f"Total_PES_Samples\n{npes}")

    logger.info(f"Completed {spin} {mol_idx_str}")


if __name__ == "__main__":
    main()
