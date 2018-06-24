"""
Testing the funnel hypothesis
"""
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import chemcoord as cc
import numpy as np
import pandas as pd

pdb = PDBFile('1ubq.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
pdb_bonds = modeller.topology.bonds()
atoms = modeller.topology.atoms()
positions = modeller.getPositions()

cc_bonds = {}
cc_positions = numpy.zeros((3, modeller.topology.getNumAtoms()))
atom_names = []

# Construct bond dictionary and positions chemcoord 
for index, atom in enumerate(atoms):
    cc_bonds[index] = set()
    pos = positions[index] / nanometer
    atom_names.append(atom.name)
    cc_positions[:, index] = pos

for bond in pdb_bonds:
    cc_bonds[bond[0].index].add(bond[1].index)

with pd.option_context('display.max_rows', None):
    print(cc_positions)
    print(atom_names)
    cc_df = pd.DataFrame({
            'atom': atom_names, 
            'x': cc_positions[0, :], 
            'y': cc_positions[1, :], 
            'z': cc_positions[2, :]
    })
    molecule = cc.Cartesian(cc_df)
    molecule.set_bonds(cc_bonds)
    zmat = molecule.get_zmat()
    molecule2 = zmat.get_cartesian()
    print(molecule2 - molecule)





