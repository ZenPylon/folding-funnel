# Exploring the Folding Funnel

(Will update later with more info description )


## TODO

1. Calculate RMSD of post-minimization (PM) protein coordinates compared with native state.
1. Write unit test for RMSD - compare biopython calculation with OpenMM calculation.
1. Calculate the RMSD between PM native state vs native state.  The RMSD should be small.
1. Setup cloud infrastructure for faster sampling.
1. Perform PM energy sampling of randomly assigned torsions.  The native state should be the smallest.
1. Plot PM angle distance vs. PM RMSD.  There should be a positive correlation.
1. Plot [RMSD between PM random samples to PM native state] vs. [energy level].
