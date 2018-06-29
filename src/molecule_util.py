from simtk.openmm.app import PDBFile, Simulation, ForceField, Modeller, PME, HBonds
from simtk.openmm import LangevinIntegrator
from simtk.openmm.vec3 import Vec3
from simtk.unit import kelvin, nanometer, picosecond, picoseconds
import numpy as np
import pandas as pd
import chemcoord as cc


class MoleculeUtil(object):
    """
    A class for managing a molecule defined by a PDB file
    """
    np.random.seed(20)

    def __init__(self, pdb_path, offset_size=2):
        # OpenMM init
        self.pdb_path = pdb_path
        self.pdb = PDBFile(self.pdb_path)
        self.forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        self.modeller = Modeller(self.pdb.topology, self.pdb.positions)
        self.modeller.addHydrogens(self.forcefield)
        self.system = self.forcefield.createSystem(
            self.modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer,
            constraints=HBonds
        )
        self.integrator = LangevinIntegrator(
            300*kelvin, 1/picosecond, 0.002*picoseconds)
        self.simulation = Simulation(
            self.modeller.topology, self.system, self.integrator)

        self.pdb_bonds = self.modeller.topology.bonds()
        self.pdb_atoms = self.modeller.topology.atoms()
        self.pdb_positions = self.modeller.getPositions()

        # Zmat and torsions
        self.cc_bonds = {}
        self.set_cc_positions(self.pdb_positions)
        self.offset_size = offset_size
        self.torsion_indices = self._get_torsion_indices()
        self.starting_torsions = np.array([
            self.zmat.loc[self.torsion_indices[:, 0], 'dihedral'],
            self.zmat.loc[self.torsion_indices[:, 1], 'dihedral']]).T
        self.seed_offsets()

    def seed_offsets(self):
        self.offsets = np.random.choice(
            [0, 0, -1, 1], self.starting_torsions.shape)

    def get_torsions(self):
        return np.array([
            self.zmat.loc[self.torsion_indices[:, 0], 'dihedral'],
            self.zmat.loc[self.torsion_indices[:, 1], 'dihedral']]).T

    def set_torsions(self, new_torsions):
        self.zmat.safe_loc[self.torsion_indices[:, 0],
                           'dihedral'] = new_torsions[:, 0]
        self.zmat.safe_loc[self.torsion_indices[:, 1],
                           'dihedral'] = new_torsions[:, 1]

    def get_offset_torsions(self, scale_factor):
        """
        Calculates and returns new torsion angles based on randomly generated
        offsets.

        Args:
            scale_factor: the relative scale of the offset relative to
                          self.offset_size

        Returns:
            The new torsion angles
        """
        total_offset = self.offset_size * scale_factor
        new_torsions = np.zeros(shape=self.starting_torsions.shape)
        new_torsions[:, 0] = self.starting_torsions[:, 0] + \
            (self.offsets[:, 0] * total_offset)
        new_torsions[:, 1] = self.starting_torsions[:, 1] + \
            (self.offsets[:, 1] * total_offset)
        return new_torsions

    def run_simulation(self):
        """
        Run a simulation to calculate the current configuration's energy level.
        Note that the atoms will likely move somewhat during the calculation,
        since energy minimization is used.

        Returns:
            A tuple of the form (potential_energy, updated_positions)
        """
        # Delete solvent that's based on previous positions
        self.modeller.deleteWater()
        cartesian = self.zmat.get_cartesian().sort_index()
        self.simulation.context.setPositions(
            [Vec3(x, y, z) for x, y, z in zip(
                cartesian['x'], cartesian['y'], cartesian['z'])]
        )
        self.modeller.addSolvent(self.forcefield, padding=1.0*nanometer)
        # self.simulation.minimizeEnergy(maxIterations=100)
        state = self.simulation.context.getState(
            getEnergy=True, getPositions=True)
        p_energy = state.getPotentialEnergy()
        positions = state.getPositions(asNumpy=True)
        print(p_energy)
        return p_energy, positions

    def _init_pdb_bonds(self):
        for index, atom in enumerate(self.pdb_atoms):
            self.cc_bonds[index] = set()

        for bond in self.pdb_bonds:
            self.cc_bonds[bond[0].index].add(bond[1].index)
            self.cc_bonds[bond[1].index].add(bond[0].index)

    def set_cc_positions(self, positions):
        """
        Calculates the zmat from an OpenMM modeller
        """
        cc_positions = np.zeros((3, self.modeller.topology.getNumAtoms()))
        atom_names = []

        # Construct bond dictionary and positions chemcoord
        for index, atom in enumerate(self.pdb_atoms):
            pos = positions[index] / nanometer
            atom_names.append(atom.name)
            cc_positions[:, index] = pos

        cc_df = pd.DataFrame({
            'atom': atom_names,
            'x': cc_positions[0, :],
            'y': cc_positions[1, :],
            'z': cc_positions[2, :]
        })

        self.cartesian = cc.Cartesian(cc_df)
        self.cartesian.set_bonds(self.cc_bonds)
        self.cartesian._give_val_sorted_bond_dict(use_lookup=True)
        self.zmat = self.cartesian.get_zmat(use_lookup=True)

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
