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
    Rotates each polypeptide backbone atom defined by a torsion angle.
    Args:
        angles (list): the angles to modify each atom.  There should be a total of
            3 * num_residues - 3 elements in the list.  The first three atoms ac
        polypeptide (list) 
    
    for i in range(0, 1):
        res = polypeptide[i]
        next_res = polypeptide[i + 1]

        try:
            a = res['N'].get_vector()
            b = res['CA'].get_vector()
            c = res['C'].get_vector()
            d = next_res['N'].get_vector()
        

        except:
            print('ERROR: missing backbone atoms\nExiting...')
            exit()