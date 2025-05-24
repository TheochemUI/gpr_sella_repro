#!/usr/bin/env python3
# TODO(rg): So much, make this a CLI

import os
import configparser

import ase.io
from ase.io import read as aseread


def mk_ef_saddle():
    _conf = configparser.ConfigParser()
    _conf.read(["config.ini"])
    _conf["Main"]["job"] = "point"
    _conf["ASE_NWCHEM"]["nwchem_path"] = os.environ["NWCHEM_COMMAND"]
    _conf["ASE_NWCHEM"]["nproc"] = "4"
    with open("config.ini", "w") as ef_conf:
        _conf.write(ef_conf)


def prep_con():
    atoms = aseread("saddle.con")
    # Centering to handle gh-188 in eOn
    atoms.set_cell([25, 25, 25])
    atoms.center()
    ase.io.write("pos.con", atoms)


if __name__ == "__main__":
    mk_ef_saddle()
    prep_con()
