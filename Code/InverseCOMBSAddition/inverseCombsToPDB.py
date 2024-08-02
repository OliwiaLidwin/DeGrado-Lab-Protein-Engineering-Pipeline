import os
import math
from pymol import cmd

# 
# Change line 46 according to the name you want your files to have.
# Replace resname 38E with the name of your ligand.
# 

def replaceResidue(file_path, target_file, output_dir, counter):
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.load(target_file, "IC_changes")
    cmd.remove("chain Y")
    cmd.select("allres", "chain X or structure")
    print(cmd.count_atoms("allres"))
    cmd.extract("binded", "allres") #Create new object.

    
    cmd.remove("structure or (resname 38E and name H1)") #No longer needed
    cmd.remove("IC_changes")
    
    cmd.select("ligand", "chain X")
    print(cmd.count_atoms("ligand"))
    
    model = cmd.get_model("ligand")
    binder_residue = model.atom[0]
    
    #Find the closest CA atom, and get it's residue
    cmd.select("CA_IC", "ligand and name CA")
    cmd.select("original_residue", "not ligand and byres name CA within 2 of CA_IC")
    
    atom_info = cmd.get_model('original_residue')
    ogResNum = atom_info.atom[0].resi
    ogResType =atom_info.atom[0].resn 
    
    cmd.remove("original_residue")
    
    cmd.alter("ligand", "chain='A'")
    cmd.alter("ligand", "segi='A'")
    cmd.alter("ligand", f"resi={ogResNum}")
    
    cmd.sort()
    cmd.rebuild()

    
    output_filename = f"Pose5_{ogResNum}_{ogResType}_to_{binder_residue.resn}_run_{counter}.pdb"
    output_path = os.path.join(output_dir, output_filename)
    cmd.save(output_path, "binded")
    return output_path