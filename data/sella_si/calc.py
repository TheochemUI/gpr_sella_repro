#!/usr/bin/env python3

import os
from ase.io import read
from ase.calculators.nwchem import NWChem
from ase.calculators.socketio import SocketIOCalculator
from sella import Sella

# Multiplicity of the molecule. 1 for singlets, 2 for doublets.
mult = 1

# Index of molecule to optimize
index = '000'

# Path to where calculations will be run
run_path = 'runs'

# Path to local scratch directory for NWChem scratch files
scratch_path = '/scratch/ehermes'

# Full path to your nwchem executable
nwchem_path = '/home/ehermes/build/nwchem/bin/LINUX64/nwchem'

# Memory for NWChem to allocate. You probably don't need to change this.
memory = '2 gb'

# Calculate the vibrational frequencies to verify final structure is a saddle point?
check_freqs = False

# Do not edit below this line

spin = 'singlets' if mult == 1 else 'doublets'

run_dir = os.path.join(run_path, spin, index)
os.makedirs(run_dir, exist_ok=True)

label = f'{spin}_{index}'
scratch_dir = os.path.join(scratch_path, label)
os.makedirs(scratch_dir, exist_ok=True)

nwchem_kwargs = dict(
    label=label,
    directory=run_dir,
    perm=scratch_dir,
    scratch=scratch_dir,
    command=f'{nwchem_path} PREFIX.nwi > PREFIX.nwo',
    memory=memory,
    scf=dict(
        nopen=mult - 1,
        thresh=1e-8,
        maxiter=200,
    ),
    basis='3-21G',
    task='optimize',
    driver=dict(
        socket=f'unix {label}',
    ),
)

if mult == 2:
    nwchem_kwargs['scf']['uhf'] = None  # switch to unrestricted calculation

nwchem = NWChem(**nwchem_kwargs)

atoms = read(f'{spin}/{index}.xyz')

with SocketIOCalculator(nwchem, unixsocket=label) as calc:
    atoms.calc = calc
    opt = Sella(
        atoms,
        internal=True,
        logfile=os.path.join(run_dir, f'{label}.log'),
        trajectory=os.path.join(run_dir, f'{label}.traj'),
    )
    opt.run(fmax=0.01, steps=200)

if check_freqs:
    nwchem_kwargs['task'] = 'frequencies'
    nwchem_kwargs['label'] = f'{label}_frequencies'
    del nwchem_kwargs['scf']['thresh']
    del nwchem_kwargs['scf']['maxiter']
    del nwchem_kwargs['driver']
    atoms.calc = NWChem(**nwchem_kwargs)
    atoms.get_potential_energy()
