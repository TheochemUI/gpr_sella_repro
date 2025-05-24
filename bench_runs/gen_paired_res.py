import ase
import ase.io as aseio
import pprint as pp
import numpy as np
import subprocess
import dataclasses
from dataclasses import dataclass, field
from pathlib import Path
from enum import Enum, StrEnum
import pandas as pd
from chemparseplot.analyze.use_ira import calculate_rmsd

gitroot = Path(
    str(
        subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True)
        .stdout.decode("utf-8")
        .strip()
    )
)

SINGLET_IDS = [f"{x:003d}" for x in range(265)]
DOUBLET_IDS = [f"{x:003d}" for x in range(235)]

class RunType(StrEnum):
    GPRD = "gprd/final_gprd_wparam"
    SELLA = "final_sella"
    ID = "idimer/final_lbfgsrot_lbfgs"


class Spin(StrEnum):
    SINGLET = "singlets"
    DOUBLET = "doublets"


class SCF_THRESH(StrEnum):
    E8M = "1e8m"
    E5M = "1e5m"


# HACK: Recall that .con files are centered in a 25 length box..
# Centering to handle gh-188 in eOn
def _set_box_center(atm: ase.Atoms):
    atm.set_cell([25, 25, 25])
    atm.center()
    return atm

def _no_ghost(atm: ase.Atoms):
    # Sella writes out X symbol'd atoms for ghost atoms
    del atm[[atom.index for atom in atm if atom.symbol=='X']]
    return atm

@dataclass
class RunPaths:
    mol_id: int
    spin: Spin
    method: RunType
    scf_thresh: SCF_THRESH = SCF_THRESH.E8M
    root_dir: Path = gitroot / "bench_runs"
    # Handled in sub-classes
    init: ase.Atoms = field(init=False)
    saddle: ase.Atoms = field(init=False)
    fdir: Path = field(init=False)


@dataclass
class eonRun(RunPaths):
    def __post_init__(self):
        assert self.method in (RunType.GPRD, RunType.ID)
        self.fdir = (
            self.root_dir
            / "eon/runs/"
            / self.scf_thresh
            / self.method
            / self.spin
            / f"{self.mol_id:003d}"
        )
        self.init = aseio.read(self.fdir / "pos.con")
        if (self.fdir / "saddle.con").exists():
            self.saddle = aseio.read(self.fdir / "saddle.con")


@dataclass
class sellaRun(RunPaths):
    def __post_init__(self):
        assert self.method == RunType.SELLA
        self.fdir = (
            self.root_dir
            / "sella/runs/"
            / self.scf_thresh
            / self.method
            / self.spin
            / f"{self.mol_id:003d}"
        )
        _traj = aseio.read(
            self.fdir / f"{self.spin}_{self.mol_id:003d}.traj", index=":"
        )
        self.init, self.saddle = _no_ghost(_set_box_center(_traj[0])), _no_ghost(_set_box_center(_traj[-1]))

all_rmsd_sella = []
for sing in SINGLET_IDS:
    sads = sellaRun(int(sing), Spin.SINGLET, RunType.SELLA).saddle
    try:
        sade = eonRun(int(sing), Spin.SINGLET, RunType.GPRD).saddle
    except Exception as _:
        continue
    all_rmsd_sella.append({'system_id': f"s_{sing}", 'sella_gprd_rmsd': calculate_rmsd(sads, sade, 10)})
for doub in DOUBLET_IDS:
    sads = sellaRun(int(doub), Spin.DOUBLET, RunType.SELLA).saddle
    try:
        sade = eonRun(int(doub), Spin.DOUBLET, RunType.GPRD).saddle
    except Exception as _:
        continue
    all_rmsd_sella.append({'system_id': f"d_{doub}", 'sella_gprd_rmsd': calculate_rmsd(sads, sade, 10)})
df = pd.DataFrame(all_rmsd_sella)

df.to_csv("sella_gprd_rmsd.csv", index=False)
print(
    f"Found {len(df)} duplicate systems. Results saved to sella_gprd_rmsd.csv"
)
