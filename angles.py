from copy import deepcopy
from Bio.PDB import *
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import Polypeptide

def rot_atom(angle: float, atoms: tuple) -> Vector:
    """
    Rotates the fourth atom "d" in a tuple of four consecutive atoms in a chain.
    Returns the new position of d.

    Algorithm taken from Practical Conversion from Torsion Space to Cartesian
    Space for In Silico Protein Synthesis
    """
    a = atoms[0].get_vector()
    b = atoms[1].get_vector()
    c = atoms[2].get_vector()
    d = atoms[3].get_vector()
    bc = c - b
    ab = b - a
    bc_normed = bc.normalized()
    
    # Start at c, and move along the bc vector (amount: length of bond bc)
    length_cd = (d - c).norm()
    d0 = bc_normed ** length_cd

    # Rotate to d1
    n = (ab ** bc_normed).normalized()
    bond_angle = calc_angle(b, c, d0)
    bond_rot = rotaxis(bond_angle, n)
    d1 = d0.left_multiply(bond_rot)

    # Rotate around bc vector based on dihedral angle
    torsion_angle = calc_dihedral(a, b, c, d)
    torsion_rot = rotaxis(torsion_angle, bc_normed)
    d2 = c + d1.left_multiply(torsion_rot)
    
    return d2

def rot_backbone(angles: list, polypeptide: Polypeptide):
    """
    Rotates each polypeptide backbone atom defined by a set of torsion angles
    Returns the newly constructed polypeptide (does not modify polypeptide param).

    Args:
        angles (list): a list of tuples containing dihedral deltas (phi, psi)
        polypeptide (list): the polypeptide to operate on
    """
    num_residues = len(polypeptide)
    if num_residues < 2:
        print('Attempting to operate on empty polypeptide.  Skipping.')
        return

    new_polypeptide = deepcopy(polypeptide)

    # This is very similar to Polypeptide.get_phi_psi_list()
    for i in range(0, num_residues):
        # TODO construct the new polypeptide
        res = new_polypeptide[i]
        n = res['N'].get_vector()
        ca = res['CA'].get_vector()
        c = res['C'].get_vector()
        
        # Phi angle
        if i > 0:
            prev_res = new_polypeptide[i - 1]
            c_prev = prev_res['C']
            new_coord = rot_atom(angles[i][0], (c_prev, n, ca, c))   
            res['C'].set_coord(new_coord)

        # (Skip the dihedral angle centered around C-N bond - assumed planar)

        # Psi angle
        if i < num_residues - 1:
            next_res = new_polypeptide[i - 1]
            n_next = next_res['C']
            new_coord = rot_atom(angles[i][1], (n, ca, c, n_next))
            next_res['N'].set_coord(new_coord)

