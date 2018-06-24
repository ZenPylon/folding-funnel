"""
Testing the funnel hypothesis
"""
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import chemcoord as cc

pdb = PDBFile('1ubq.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
bonds = modeller.topology.bonds()

for bond in bonds:
    print(bond)