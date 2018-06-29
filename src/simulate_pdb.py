from simtk.openmm.app import *

from simtk.unit import *
from sys import stdout
import json


state = simulation.context.getState(getEnergy=True)
print(state.getPotentialEnergy())
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
