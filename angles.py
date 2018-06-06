from Bio.PDB import *

def rot_atom(angle: float, atoms: tuple):
    """
    Rotates the fourth atom "d" in a tuple of four consecutive atoms in a chain.
    Returns the new position of d.

    Algorithm taken from Practical Conversion from Torsion Space to Cartesian
    Space for In Silico Protein Synthesis
    """
    a = atoms[0], b = atoms[1], c = atoms[2], d = atoms[3]
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

def rot_backbone(angles: list, polypeptide: polypeptide):
    """
    Rotates each polypeptide backbone atom defined by a set of torsion angles
    Returns the newly constructed polypeptide (does not modify polypeptide param).

    Args:
        angles (list): the angles to modify each atom.  There should be a total of num_residues - 2 angles
        polypeptide (list): the polypeptide to operate on
    """
    # TODO get the first n, ca, c, n bond here

    for i in range(1, len(polypeptide - 1)):
        # TODO construct the new polypeptide
        res = polypeptide[i]
        next_res = polypeptide[i + 1]
        a = res['N'].get_vector()
        b = res['CA'].get_vector(),
        c = res['C'].get_vector()
        d = next_res['N'].get_vector()
        e = next_res['CA'].get_vector()
        f = next_res['C'].get_vector()

        rot_d = rot_atom(angle, (a, b, c, d))
        rot_e = rot_atom(angle, (b, c, d, e))
        rot_f = rot_atom(angle, (c, d, e, f))
    

    # TODO get the last c, n, ca, c bond here
