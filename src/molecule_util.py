from simtk.openmm.app import PDBFile, ForceField, Modeller, PME, HBonds
from simtk.unit import kelvin, nanometer, picosecond, picoseconds
import numpy as np
import pandas as pd
import chemcoord as cc


class MoleculeUtil(object):
    """
    A class for managing a molecule defined by a PDB file
    """
    np.random.seed(20)

    def __init__(self, pdb_path, offset_size=4):
        self.offset_size = offset_size
        self.modeller = self._get_modeller(pdb_path)
        self.zmat = self._get_zmat()
        self.torsion_indices = self._get_torsion_indices()
        self.starting_torsions = np.array([
                self.zmat.loc[self.torsion_indices[:, 0], 'dihedral'],
                self.zmat.loc[self.torsion_indices[:, 1], 'dihedral']]).T
        print(self.starting_torsions)

    def get_new_torsions(self, scale_factor):
        """
        Calculates and returns new torsion angles based on randomly generated
        offsets.

        Args:
            scale_factor: the relative scale of the offset relative to 
                          self.offset_size

        Returns:
            The new torsion angles
        """
        offsets = np.random.choice([0, 0, -1, 1], self.starting_torsions.shape)
        total_offset = self.offset_size * scale_factor
        new_torsions[:, 0] = starting_torsions[:, 0] + (offsets[:, 0] * total_offset)
        new_torsions[:, 1] = starting_torsions[:, 1] + (offsets[:, 1] * total_offset)
        return new_torsions
        
    def _get_modeller(self):
        pdb = PDBFile(self.pdb_file)
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield)
        return modeller

    def _get_zmat(self):
        """
        Calculates the zmat from an OpenMM modeller
        """
        # Create new document with field pdb_name, doc_id as random string
        pdb_bonds = self.modeller.topology.bonds()
        atoms = self.modeller.topology.atoms()
        positions = self.modeller.getPositions()

        cc_bonds = {}
        cc_positions = np.zeros((3, modeller.topology.getNumAtoms()))
        atom_names = []

        # Construct bond dictionary and positions chemcoord
        for index, atom in enumerate(atoms):
            cc_bonds[index] = set()
            pos = positions[index] / nanometer
            atom_names.append(atom.name)
            cc_positions[:, index] = pos

        for bond in pdb_bonds:
            cc_bonds[bond[0].index].add(bond[1].index)
            cc_bonds[bond[1].index].add(bond[0].index)

        cc_df = pd.DataFrame({
            'atom': atom_names,
            'x': cc_positions[0, :],
            'y': cc_positions[1, :],
            'z': cc_positions[2, :]
        })

        molecule = cc.Cartesian(cc_df)
        molecule.set_bonds(cc_bonds)
        molecule._give_val_sorted_bond_dict(use_lookup=True)
        zmat = molecule.get_zmat(use_lookup=True)
        return zmat

    def _get_torsion_indices(self):
        """
        Calculates indices into the zmatrix which correspond to phi
        and psi angles.

        Args:
            zmat: the zmatrix specifying the molecule
        Returns:
            a numpy.array, with first column as phi_indices, second column 
            as psi_indices
        """
        phi_indices = []
        psi_indices = []

        for i in range(len(self.zmat.index)):
            b_index = self.zmat.loc[i, 'b']
            a_index = self.zmat.loc[i, 'a']
            d_index = self.zmat.loc[i, 'd']

            # If this molecule references a magic string (origin, e_x, e_y, e_z, etc)
            if isinstance(b_index, str) or isinstance(a_index, str) or isinstance(d_index, str):
                continue

            # Psi angles
            if (self.zmat.loc[i, 'atom'] == 'N') & \
                    (self.zmat.loc[b_index, 'atom'] == 'CA') & \
                    (self.zmat.loc[a_index, 'atom'] == 'C') & \
                    (self.zmat.loc[d_index, 'atom'] == 'N'):
                psi_indices.append(i)

            elif (self.zmat.loc[i, 'atom'] == 'N') & \
                    (self.zmat.loc[b_index, 'atom'] == 'C') & \
                    (self.zmat.loc[a_index, 'atom'] == 'CA') & \
                    (self.zmat.loc[d_index, 'atom'] == 'N'):
                psi_indices.append(i)

            elif (self.zmat.loc[i, 'atom'] == 'C') & \
                    (self.zmat.loc[b_index, 'atom'] == 'N') & \
                    (self.zmat.loc[a_index, 'atom'] == 'CA') & \
                    (self.zmat.loc[d_index, 'atom'] == 'C'):
                phi_indices.append(i)

            elif (self.zmat.loc[i, 'atom'] == 'C') & \
                    (self.zmat.loc[b_index, 'atom'] == 'CA') & \
                    (self.zmat.loc[a_index, 'atom'] == 'N') & \
                    (self.zmat.loc[d_index, 'atom'] == 'C'):
                phi_indices.append(i)

        return np.array([phi_indices, psi_indices]).T
