import util

# Don't worry about hydrogen atoms - we'll add them to the PDB model later
structure, residue_list, polypeptide = util.main_load('ubiq', '1ubq.pdb')
# torsion_angles = polypeptide.get_phi_psi_list()
# backbone_vectors = []
# res_count = len(polypeptide)



# TODO - refactor logic into function, then apply function to every string of 4 atoms
# TODO - test by comparing coordinates of original with constructed coordinates
#        (without modifying the coordinates)
# TODO - test by changing angles more and more, and plot vs. RMSD distance
#        Expectation is that it RMSD should increase as angle distance increases



# print(backbone_vectors)


# Construct backbone from angles and distances

# 1. Get torsion angles and bond distances for backbone molecules
# 1. Get vectors for sidechain atoms relative to backbone atom
# 1. Iterate through each residue and construct each backbone atom 
#    relative to transformed bond angle (use spherical coordinates)


# THEN
# 1. Save PDB with minor alterations
# 1. Calculate RMSD.  Should be small
# 1. Plot superimposed or side-by-side structures.  Get visual confirmation
# of similarity

# THEN:
# 1. Make increasingly large deviations from loaded structure
# 1. Calculate RMSD
# 1. Plot Angle distance vs RMSD.  Should be positive correlation



# TODO
# Try to use ChemCoords to convert between coordinates and bonds 
# (given by pdb positions and topology).  Verify that
# 1. Transforming between internal and cartesian stays the same when 
#    no modifications are made
# 2. Modifying the angle of a residue halfway through the chain only modifies 
#    coordinates on that residue or after (no residues before).

# If ChemCoords is unworkable, just construct the backbone iteratively through 
# local spherical coordinate --> cartersian coordinate geometry, and keep the side
# chains in the same relative position (relative to the alpha carbon).

# TODO: The Algorithm:
# 
# 1. Specify all phi / psi angle deltas (magnitude M) relative to native structure
#    and get absolute phi / psi coordinates (calculate native angles + deltas).
# 2. Convert from angle-space to cartesian coordinates.
# 3. Pass coordinates to OpenMM (how?).
    # Use topology.addChain, topology.addResidue(), etc based on constructed structure
# 4. Add solvent to model with OpenMM.
# 5. Run energy minimization on model.
# 6. Pass final energy value and coordinates to BioPython, and reset atom coords.
# 7. Recalculate phi / psi angles from coordinates.
# 8. Calculate distance in angle-space from starting point (in this case,
#    native structure) and calculate delta in energy level.
# 9. Output distance and energy-delta to file (i.e track progress).
# 10. Repeat steps 1-9, increasing delta magnitude M each time.
