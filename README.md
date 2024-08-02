Pipeline for protein-ligand binding from COMBS output to MD Simulation Analysis.
A) Planarity 
      - Finding how planar the ligand is with the initial interacting residues.
      - Creating a bar plot highlighting all the outputs that fit the criteria.

B) inverseCOMBSAddition 
      - Addition of inverseCOMBS output to COMBS output.

C) inverseCOMBStoScripts
      - Generate the rest of the protein sequence using binding site context and with ligand context (LigandMPNN) or without consideration for the ligand context (Chroma, RFDiffusion)

D) RMSD_ESMFold
      - Find side chain and backbone rmsd. 
      - Add the ligand back into the ESMFold output
      - Create bar plots for best rmsd valued proteins.

E) CPPTRAJ (needs to be cited)
      - Analyze MD simulation output by dihedral angle, backbone RMSD, and hydrogen bonding.
