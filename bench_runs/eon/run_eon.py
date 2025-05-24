#!/usr/bin/env python3

import os
from io import StringIO
import shutil
from pathlib import Path
from zipfile import ZipFile
import configparser
import logging

import ase.io
from ase.io import read as aseread

from rgpycrumbs._aux import getstrform, get_gitroot, switchdir

import eon.akmc
from eon.explorer import MinModeExplorer
from eon.config import ConfigClass as eonConf
import eon.fileio as eio

import click

# XXX: eon.config.config is a leaky global state.. why.


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
    ase.io.write(getstrform(run_dir / "pos.con"), atm)

def prepare_con_file(run_dir, spin, mol_idx_str, gitroot):
    """Prepares the pos.con file."""
    with ZipFile(gitroot / "data" / "sella_si_data.zip", "r") as zdat:
        with zdat.open(f"sella_si/{spin}/{mol_idx_str}.xyz", "r") as atmdat:
            atoms = aseread(StringIO(atmdat.read().decode()), format="xyz")
            # Centering to handle gh-188 in eOn
            atoms.set_cell([25, 25, 25])
            atoms.center()
            ase.io.write(getstrform(run_dir / "pos.con"), atoms)


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


def prepare_config_ini(run_dir, config_type, nwchem_path, mult, scfthresh):
    """Prepares the config.ini file."""
    _conf = configparser.ConfigParser()
    config_files = ["base_config.ini"]

    if config_type == "gprd":
        config_files.append("gpr_mixin.ini")
        _conf.read(config_files)
        _conf["Saddle Search"]["min_mode_method"] = "gprdimer"
    elif config_type == "idimer":
        config_files.append("dimer_mixin.ini")
        _conf.read(config_files)
        _conf["Saddle Search"]["min_mode_method"] = "dimer"
    elif config_type == "min":
        _conf.read(config_files)
        _conf["Main"]["job"] = "minimization"
    else:
        raise ValueError("Invalid config_type specified.")

    _conf["ASE_NWCHEM"]["nwchem_path"] = nwchem_path
    _conf["ASE_NWCHEM"]["multiplicity"] = str(mult)
    _conf["ASE_NWCHEM"]["nproc"] = os.environ["NWCHEM_MPI"]
    _conf["ASE_NWCHEM"]["scf_thresh"] = str(convert_scfthresh(scfthresh))

    with open(run_dir / "config.ini", "w") as akmc_conf:
        _conf.write(akmc_conf)


def setup_and_run_initial_displacement(run_dir, logger):
    """Sets up and runs the initial displacement."""
    with switchdir(run_dir):
        econf = eonConf()
        econf.init("config.ini")

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
    default="gprd_norm_fixed",
    help="Run directory name within runs/{scfthresh}/",
)
@click.option(
    "--config-type",
    type=click.Choice(["gprd", "idimer", "min"]),
    default="gprd",
    help="Type of config to use: gpr_mixin.ini or dimer_mixin.ini",
)
@click.option(
    "--min_dir",
    type=click.Path(
        exists=False, file_okay=False, dir_okay=True, readable=True, path_type=Path
    ),
    default=None,
    help="Directory containing a min.con file to use instead of generating pos.con.",
)
def main(mol_idx, spin, scfthresh, rdir, config_type, min_dir):
    """GPRD benchmarks on the Sella test set.

    MOL_IDX is the index of the molecule to optimize (e.g., 011).
    """
    logger = setup_logger()

    gitroot = get_gitroot()
    mult = 1 if spin == "singlets" else 2
    mol_idx_str = f"{mol_idx:03d}"
    run_path = Path.cwd() / "runs" / scfthresh / config_type / rdir
    nwchem_path = os.environ["NWCHEM_COMMAND"]

    run_dir = run_path / spin / mol_idx_str
    run_dir.mkdir(parents=True, exist_ok=True)

    if min_dir:
        prepare_min_con(
            run_dir,
            spin,
            mol_idx_str,
            min_dir,
        )
    else:
        prepare_con_file(run_dir, spin, mol_idx_str, gitroot)
    prepare_config_ini(run_dir, config_type, nwchem_path, mult, scfthresh)
    setup_and_run_initial_displacement(run_dir, logger)

    logger.info(f"Preparing {spin} {mol_idx_str} with config type: {config_type}")


if __name__ == "__main__":
    main()
