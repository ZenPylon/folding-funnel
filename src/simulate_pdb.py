from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import json

pdb = PDBFile('1ubq.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
bonds = modeller.topology.bonds()
# simulation = Simulation(modeller.topology, system, integrator)
# simulation.context.setPositions(modeller.positions)


# Get all bonds between atoms in the same residue and the polypeptide bond
# between alpha carbons and the next nitrogen
for bond in bonds:
    # inner_bond = bond[0].residue.index == bond[1].residue.index
    # peptide_bond = 
    # if inner_bond :
    print(bond)



# simulation.minimizeEnergy(maxIterations=100)
# state = simulation.context.getState(getEnergy=True)
# print(state.getPotentialEnergy())
# simulation.reporters.append(PDBReporter('output.pdb', 1000))
# simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
# simulation.step(10000)
