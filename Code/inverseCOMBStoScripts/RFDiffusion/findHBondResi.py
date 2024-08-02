from pymol import cmd

# 
# Find all of the hydrogen bonding residues based on a certain distance.
# Change line 11 accordingly.
# 
def selectE(file_path, ligandName):
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.select("ligand", f"resname {ligandName}")
    cmd.select("interacting_residues", "not ligand and chain A and byres (elem O or elem N within 3.4 of (ligand and elem O or ligand and elem N))")    
    NumberofAtoms = cmd.count_atoms('interacting_residues')
    if NumberofAtoms>1:
        print(f"Number of atoms near ligand: {NumberofAtoms}")
    atom_info = cmd.get_model("interacting_residues")
    residues = []
    for atom in atom_info.atom:
        residues.append(atom.resi)
    residues = list(set(residues))
    num_residues = len(residues)
    print(f"Number of residues near {ligandName}: {num_residues}")
    return residues

def numberOfResidues(file_path, ligandName=""):
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.select("ligand", f"resname {ligandName}")
    cmd.select("protein", "not ligand")
    atom_info = cmd.get_model("protein")
    residues = []
    for atom in atom_info.atom:
        residues.append(atom.resi)
    residues = list(set(residues))
    print(len(residues))
    return len(residues)

def hasLigand(file_path, ligandName):
    print(ligandName)
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.select("ligand", f"resname {ligandName}")
    print(cmd.count_atoms("ligand"))
    if cmd.count_atoms("ligand") == 0:
        return False
    return True