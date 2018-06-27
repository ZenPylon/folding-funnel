# Exploring the Folding Funnel

(Will update later with more info description )


## Algorithm

### Messages and Data
The response data in `init_dispatch` contains:
- The starting torsion angles (phi, psi) list

The response data in `job_dispatch` contains:
- Modified angles (phi, psi) list

### Process
1. Start broker process
1. Broker process downloads PDB from cloud storage
1. Broker initializes zmatrix from PDB topology
1. Broker process listens for worker `init_request` 
1. Broker process listens for worker `job_request`
1. Start worker process(es)
1. Worker downloads PDB file from cloud storage
1. Worker initializes OpenMM simulation from PDB file data
1. Worker initializes biopython chain from PDB file
1. Worker initializes zmatrix from PDB topology
1. Worker process listens for broker `init_dispatch`
1. Worker process listens for broker `job_dispatch`
1. Worker process sends `init_request` message
1. Broker process responds with `init_dispatch` message/data
1. Worker process sends `job_request` message
1. Broker process sends `job_dispatch` message/data
1. Worker converts modified angles back into cartesian coordinates using chemcoord
1. Worker constructs OpenMM Topology with updated cartesian coordinates
1. Worker runs simulation to minimize energy level
1. Worker records final energy level and final coordinates upon simulation completion
1. Worker constructs biopython chain with updated coordinates
1. Worker calculates RMSD between new chain and original chain using biopython
1. Worker recalculates torsions with chemcoord based on minimized simulated/structure
1. Worker records distance in angle space
1. Worker saves distance in angle space, RMSD, and energy level to firebase
1. Worker requests a new job with `job_request` message


