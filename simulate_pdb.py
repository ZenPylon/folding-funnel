from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import json

pdb = PDBFile('1ubq.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
# modeller.addSolvent(forcefield, padding=1.0*nanometers)

system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
print(pdb.topology)

bonds = pdb.topology.bonds()

for bond in bonds:
    print(bond)
# simulation.minimizeEnergy(maxIterations=100)
# state = simulation.context.getState(getEnergy=True)
# print(state.getPotentialEnergy())
# simulation.reporters.append(PDBReporter('output.pdb', 1000))
# simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
# simulation.step(10000)
