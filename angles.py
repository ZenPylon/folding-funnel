from copy import deepcopy
from math import pi
import numpy as np
from Bio.PDB import *
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import Polypeptide


def rot_atom(torsion_angle: float, atoms: tuple) -> Vector:
    """
    Rotates the fourth atom "d" in a tuple of four consecutive atoms in a chain by
    setting to the current torsion angle to `torsion_angle`
    Returns the new position of d as a Vector.

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

    # Start at c, and move along the bc vector (amount: length of bond cd)
    length_cd = (d - c).norm()
    d0 = bc_normed ** length_cd

    # Rotate to d1
    n = (ab ** bc_normed).normalized()

    # TODO - Why is this "pi - " fudge factor needed?  What's going wrong?
    bond_angle = pi - calc_angle(b, c, d)
    bond_rot = rotaxis(bond_angle, n)
    d1 = d0.left_multiply(bond_rot)

    # Rotate around bc vector based on dihedral angle
    # torsion_angle = calc_dihedral(a, b, c, d)
    torsion_rot = rotaxis(torsion_angle, bc_normed)
    d2 = d1.left_multiply(torsion_rot)
    # print('Vectors:')
    # print(a)
    # print(b)
    # print(c)
    # print(d)
    # print(f'bc {bc}')
    # print(f'bc normed {bc_normed}')
    # print(f'ab {ab}')
    # print(f'length cd {length_cd}')
    # print(f'n {n}')
    # print(f'd0 {d0}')
    # print(f'bond angle {bond_angle}')
    # print(f'bond rot {bond_rot}')
    # print(f'd1 {d1}')
    # print(f'torsion angle {torsion_angle}')
    # print(f'torsion rot {torsion_rot}')
    # print(f'd2 {d2}')
    return c + d2


def rot_backbone(angles: list, polypeptide: Polypeptide):
    """
    Rotates each polypeptide backbone atom defined by a set of torsion angles
    Returns the newly constructed polypeptide (does not modify polypeptide param).

    Args:
        angles (list): a list of 3-tuples containing dihedral deltas (phi, psi, omega)
        polypeptide (list): the polypeptide to operate on
    """
    num_residues = len(polypeptide)
    if num_residues < 2:
        print('Attempting to operate on empty polypeptide.  Skipping.')
        return

    new_polypeptide = deepcopy(polypeptide)

    for i in range(0, num_residues):
        # TODO construct the new polypeptide
        res = new_polypeptide[i]
        n = res['N']
        ca = res['CA']
        c = res['C']

        if i > 0:
            prev_res = new_polypeptide[i - 1]
            c_prev = prev_res['C']
            ca_prev = prev_res['CA']

            # Omega angle - set this first, since psi1 is followed by omega1
            # which is followed by phi1, followed by psi2, etc.
            new_coord_ca = rot_atom(angles[i][2], (ca_prev, c_prev, n, ca))
            res['CA'].set_coord(np.array(new_coord_ca.get_array()))
            
            # Phi angle
            new_coord_c = rot_atom(angles[i][0], (c_prev, n, ca, c))
            res['C'].set_coord(np.array(new_coord_c.get_array()))

        # Psi angle
        if i < num_residues - 1:
            next_res = new_polypeptide[i + 1]
            n_next = next_res['N']
            new_coord_n = rot_atom(angles[i][1], (n, ca, c, n_next))
            next_res['N'].set_coord(np.array(new_coord_n.get_array()))


    return new_polypeptide


# Basically a copy of get_phi_psi_list, but gets omegas too
def get_all_backbone_torsions(polypeptide):
    """Return the list of phi/omega/psi dihedral angles."""
    ppl = []
    lng = len(polypeptide)
    for i in range(0, lng):
        res = polypeptide[i]
        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()
        except Exception:
            # Some atoms are missing
            # Phi/Psi cannot be calculated for this residue
            ppl.append((None, None))
            res.xtra["PHI"] = None
            res.xtra["PSI"] = None
            continue
        # Phi
        if i > 0:
            rp = polypeptide[i - 1]
            try:
                cp = rp['C'].get_vector()
                cap = rp['CA'].get_vector()
                phi = calc_dihedral(cp, n, ca, c)
                omega = calc_dihedral(cap, cp, n, ca)
            except Exception:
                omega = None
                phi = None
        else:
            # No phi for residue 0!
            omega = None
            phi = None
        # Psi
        if i < (lng - 1):
            rn = polypeptide[i + 1]
            try:
                nn = rn['N'].get_vector()
                psi = calc_dihedral(n, ca, c, nn)
            except Exception:
                psi = None
        else:
            # No psi for last residue!
            psi = None
        ppl.append((phi, psi, omega))
        # Add Phi/Psi to xtra dict of residue
        res.xtra["PHI"] = phi
        res.xtra["OMEGA"] = omega
        res.xtra["PSI"] = psi
    return ppl
