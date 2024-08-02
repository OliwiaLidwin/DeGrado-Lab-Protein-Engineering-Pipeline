import math
import sys
from pymol import cmd

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

def selectE(file_path, ligandName):
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.select("ligand", f"resname {ligandName}")
    cmd.select("interacting_residues", "not ligand and chain A and byres ((elem O and not ligand) within 4.4 of ((ligand and name O2) or (ligand and name N1) or (ligand and name O1)))")
    NumberofAtoms = cmd.count_atoms('interacting_residues')
    if NumberofAtoms>1:
        print(f"Number of atoms near ligand: {NumberofAtoms}")
    atom_info = cmd.get_model("interacting_residues")
    residues = []
    for atom in atom_info.atom:
        residues.append(atom.resi)
    residues = list(set(residues))
    residues = [int(num_str) for num_str in residues]
    residues = sorted(residues)
    num_residues = len(residues)
    print(f"Number of residues near {ligandName}: {num_residues}")
    return residues

def distance(reference_coords, esm_coords):
    distance = 0
    for r_atom,e_atom in zip(reference_coords,esm_coords):
        print(r_atom)     
        print(e_atom)      
        distance = 0      
        if r_atom[0] != e_atom[0]:
            print("ORDERING IS WRONG")
            sys.exit(1)  
        refernce_coordinates = r_atom[1]
        r_x, r_y, r_z = refernce_coordinates
        refernce_coordinates = e_atom[1]
        e_x, e_y, e_z = refernce_coordinates

        distance_cal = (math.pow((e_x-r_x),2) + math.pow((e_y-r_y),2) + math.pow((e_z-r_z),2))
        distance = distance + distance_cal
    return distance
    
def get_residue_coordinates(selection):
    atom_data = []
    for atom in selection.atom:
        atom_data.append((atom.name, (atom.coord[0], atom.coord[1], atom.coord[2])))
    sorted_atom_data = sorted(atom_data, key=lambda x: x[0])
    print(sorted_atom_data)
    return sorted_atom_data

def get_residue_resn(selection_name):
    residue_resn = None
    cmd.iterate_state(1, selection_name, 'residue_resn = resn', space=locals())
    return residue_resn

def get_residue_names(selection):
    resn_list = []
    cmd.iterate(selection, 'resn_list.append(resn)', space={'resn_list': resn_list})
    return list(set(resn_list))



def RMSD_Calc(reference_pdb, esm_pdb):
    print(reference_pdb)
    print(esm_pdb)
    rmsd_sc = 0
    rmsd_bb = 0
    residueTot = numberOfResidues(reference_pdb, '38E')
    interactingResidues = selectE(reference_pdb, '38E')
    cmd.reinitialize()

    cmd.load(reference_pdb, "referenceLigand")
    print(cmd.count_atoms("referenceLigand"))
    cmd.load(esm_pdb, "esmfold")
    print(cmd.count_atoms("esmfold"))

    print(cmd.get_object_list())
    cmd.create("ref", "referenceLigand and not resname 38E")
    print(cmd.count_atoms("ref"))
    result = cmd.cealign("esmfold", "ref")
    print(cmd.count_atoms("ref"))
    print(cmd.count_atoms("esmfold"))

    print(result['RMSD'])
    rmsd_bb = result['RMSD']
    
    # print(cmd.get_object_list())
    # cmd.remove("referenceLigand")
    # print(cmd.get_object_list())
    # print(cmd.count_atoms("ref"))
    # print(interactingResidues)
    # for residueNum in range(1, residueTot + 1):
    #     if residueNum in interactingResidues:
    #         print(cmd.get_object_list())
    #         print(residueNum)
    #         cmd.select(f"reference_resi1", f"ref and resi {residueNum} and not bb. and not elem H")
    #         print(cmd.count_atoms(f"reference_resi1"))

    #         cmd.select(f"esm_resi1", f"esmfold and resi {residueNum} and not bb. and not elem H")
    #         print(cmd.count_atoms(f"esm_resi1"))
            
            
    #         resn_ref = get_residue_names("reference_resi1")
    #         resn_esm = get_residue_names("esm_resi1")

    #         # Print residue names for debugging
    #         print(f"Reference resn: {resn_ref}")
    #         print(f"Esmfold resn: {resn_esm}")
                
        
    #         reference_model1 = cmd.get_model(f"reference_resi1")
    #         reference_coords1 = get_residue_coordinates(reference_model1)
    #         esm_model1 = cmd.get_model(f"esm_resi1")
    #         esm_coords1 = get_residue_coordinates(esm_model1)
    #         rmsd_sc = rmsd_sc + distance(reference_coords1, esm_coords1)
    #         print(rmsd_sc)
    # rmsd_sc = math.pow(((1/len(interactingResidues))*rmsd_sc), 0.5)
    return rmsd_sc, rmsd_bb
