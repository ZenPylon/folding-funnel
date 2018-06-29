from math import pi
from simtk.openmm.app import PDBFile, ForceField, Modeller, PME, HBonds
from simtk.unit import kelvin, nanometer, picosecond, picoseconds
from molecule_util import get_modeller, calculate_zmat
import chemcoord as cc
import numpy as np
import pandas as pd
import util

def run_sequence(zmat: cc.Zmat, starting_torsions, torsion_indices, offsets, offset_size=4, num_configs=3):
    """
    Args:
        pdb_name: name field in firebase
        pdb_file: path to pdb file
        offset_size: amount to deviate on each configuration
        num_configs: number of configurations to try
    """

    # Modify (on average) half the torsion angles
    for i in range(num_configs):
        zmat.safe_loc[torsion_indices[:, 0], 'dihedral'] = \
            starting_torsions[:, 0] + (offsets[:, 0] * i * offset_size)
        zmat.safe_loc[torsion_indices[:, 1], 'dihedral'] = \
            starting_torsions[:, 1] + (offsets[:, 1] * offset_size * i)
        print(zmat.loc[torsion_indices[:, 0], 'dihedral'])
        print(zmat.loc[torsion_indices[:, 1], 'dihedral'])


init_model()

