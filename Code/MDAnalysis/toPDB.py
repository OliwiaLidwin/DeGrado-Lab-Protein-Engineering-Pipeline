import MDAnalysis as mda

# Paths to your files
nc_file = 'C:/Users/oliwi/combs/Code/MDAnalysis/inputDir/Pose31_N2C8O2_2_modified_structure_resi_109_GLY_to_GLN_run_14_4/step5.nc'
parm7_file = 'C:/Users/oliwi/combs/Code/MDAnalysis/inputDir/Pose31_N2C8O2_2_modified_structure_resi_109_GLY_to_GLN_run_14_4/Pose31_N2C8O2_2_modified_structure_resi_109_GLY_to_GLN_run_14_4_0.12385695877038075_0.39356860440650354_renumbered.parm7'
output_pdb_file = 'C:/Users/oliwi/combs/Code/MDAnalysis/outputDir/run_14_4.pdb'

# Load the trajectory and topology
u = mda.Universe(parm7_file, nc_file)

# Write the trajectory to PDB format
with mda.Writer(output_pdb_file, bonds=None, n_atoms=u.atoms.n_atoms) as pdb:
    for ts in u.trajectory:
        pdb.write(u.atoms)

print(f"PDB file '{output_pdb_file}' has been created.")