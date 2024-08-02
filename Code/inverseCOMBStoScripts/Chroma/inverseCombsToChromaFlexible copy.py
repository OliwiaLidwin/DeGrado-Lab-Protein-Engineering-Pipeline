#FLEXIBLE
import os
import sys
import torch
import chroma
from pymol import cmd
from chroma import Chroma, Protein, api, conditioners
from chroma.layers.structure.conditioners import SubsequenceConditioner
from chroma.models import graph_classifier, procap, Chroma

def numberOfResidues(file_path, ligandName):
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
    cmd.select("interacting_residues", "not ligand and chain A and byres (elem O or elem N within 3.4 of (ligand and elem O or ligand and elem N))")
    NumberofAtoms = cmd.count_atoms('interacting_residues')
    if NumberofAtoms > 1:
        print(f"Number of atoms near ligand: {NumberofAtoms}")
    atom_info = cmd.get_model("interacting_residues")
    residues = []
    for atom in atom_info.atom:
        residues.append(int(atom.resi))  # Ensure residues are integers
    residues = list(set(residues))
    num_residues = len(residues)
    print(f"Number of residues near {ligandName}: {num_residues}")
    return residues

def ChromaRun(pwl_path, p_path, p_output_dir, numOutputs):
    length = numberOfResidues(pwl_path, '38E')
    mask = torch.zeros(1, length, device='cuda')
    fixed_residue = selectE(pwl_path, '38E')
    mask = torch.ones(1, length, device='cuda')
    for i in range(10):
        chroma_instance = Chroma()
        noiser=chroma.layers.structure.diffusion.DiffusionChainCov()
        mask[0, torch.tensor(fixed_residue) - 1] = 0
        protein = Protein.from_PDB(f'{p_path}', device='cuda')
        X,C,S=protein.to_XCS()
        X_noised=noiser.forward(X, C, 0.5).detach()
        noised_protein = protein.from_XCS(X_noised, C, S)
        conditioner = SubsequenceConditioner(chroma_instance.design_network, protein, mask)
        protein1=chroma_instance.sample(protein_init=noised_protein,conditioner=conditioner,tspan=[0.5, 0.001],steps=250,design_selection=mask)        
        elbo = chroma_instance.score(protein1)['elbo'].score
        output_file = os.path.join(p_output_dir, f"{os.path.basename(pwl_path)}_{i}_{elbo}.pdb")
        protein1.to(output_file)

def main():
    if len(sys.argv) < 4:
        print("Usage: python <script_name> <protein_PDB_L_dir> <protein_PDB_Dir> <output_dir_protein_PDBS> <optional_additional_arg>")
        return  # Exit the function if there are not enough arguments

    protein_L_dir = sys.argv[1]
    protein_pdb_dir = sys.argv[2]
    protein_output_dir = sys.argv[3]

    if not os.path.isdir(protein_L_dir):
        print(f"Input directory for protein PDBs '{protein_L_dir}' is not a directory.")
        return  # Exit the function if the directory is not valid

    if not os.path.isdir(protein_pdb_dir):
        print(f"Input directory for protein PDBs without ligand '{protein_pdb_dir}' is not a directory.")
        return  # Exit the function if the directory is not valid

    print("protein_pdb_dir is a valid directory")
    print("protein_pdb_dir without ligand is a valid directory")

    # Ensure output directory exists
    os.makedirs(protein_output_dir, exist_ok=True)
    api.register_key('48f0eeace33f42abaf846ccdc5f6dd30')
    
    files1 = os.listdir(protein_L_dir)
    files2 = os.listdir(protein_pdb_dir)
    files1.sort()
    files2.sort()

    for file1, file2 in zip(files1, files2):
        if not file1.endswith('.pdb'):
            print(f"Skipping non-PDB file: {os.path.basename(file1)}")
            continue
        
        file1_path = os.path.join(protein_L_dir, file1)
        file2_path = os.path.join(protein_pdb_dir, file2)
        file_output_dir = os.path.join(protein_output_dir, os.path.basename(file1))
        os.makedirs(file_output_dir, exist_ok=True)
        
        if len(sys.argv) > 4:
            additional_arg = sys.argv[4]
            ChromaRun(file1_path, file2_path, file_output_dir, additional_arg)
        else:
            ChromaRun(file1_path, file2_path, file_output_dir, None)

if __name__ == "__main__":
    main()