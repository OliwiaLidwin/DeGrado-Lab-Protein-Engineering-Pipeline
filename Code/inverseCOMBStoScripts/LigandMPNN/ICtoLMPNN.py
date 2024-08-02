#Writes and executes the scripts for ligandMPNN
import sys
import os
from pymol import cmd

# 
# Generates Scripts to run LigandMPNN with.
# This usage is for windows users as LigandMPNN needs to be run on a MacOS/Linux supported device. 
# Ubuntu was used. 
# Users should redefine ScriptProcessing accordingly.
# 

def selectE(file_path, ligandName):
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.select("ligand", f"resname {ligandName}")
    cmd.select("interacting_residues", "not ligand and chain A and byres ((elem O and not (ligand and (name O2 or name O1))) or (elem N and not(ligand and name N2)) within 4.4 of ((ligand and name O2) or (ligand and name N1) or (ligand and name O1)))")
    
    cmd.select("interacting_residues_amine", "not ligand and chain A and byres ((elem O and not (ligand and (name O2 or name O1))) or (elem N and not (ligand and name N2)) within 4.4 of (ligand and name N1))")
    NumberofAtoms = cmd.count_atoms('interacting_residues')
    NumberofAtoms2 = cmd.count_atoms('interacting_residues_amine')
    if NumberofAtoms > 1 and NumberofAtoms2 > 1:
        print(f"Number of atoms near ligand: {NumberofAtoms}")
        print(f"Number of atoms near ligand: {NumberofAtoms2}")
    atom_info = cmd.get_model("interacting_residues")
    atom_info2 = cmd.get_model("interacting_residues_amine")
    residues = []
    residues2 = []
    for atom in atom_info.atom:
        residues.append(int(atom.resi))  # Ensure residues are integers
    for atom2 in atom_info2.atom:
        residues2.append(int(atom2.resi))  # Ensure residues are integers
    
    residues = list(set(residues))
    residues2 = list(set(residues2))    
    residues = sorted(residues)
    print(residues)
    print(residues2)

    num_residues = len(residues) - len(residues2)
    print(f"Number of residues near {ligandName}: {num_residues}")
    residues = 'A' + ' A'.join(map(str, residues))
    return residues, num_residues

def ScriptProcessing(protein_path, output_directory, residues):
    input = "/mnt/c/Users/oliwi/Combs/730"
    protein_path = os.path.join(input, os.path.basename(protein_path))
    output = "/mnt/c/Users/oliwi/combs/Code/LigandMPNN/730"
    output_directory = os.path.join(output, os.path.basename(protein_path))
    script = f"""
    
python3 /mnt/c/Users/oliwi/LigandMPNN/run.py --model_type "ligand_mpnn" --seed 111 --pdb_path "{protein_path}" --out_folder "{output_directory}" --ligand_mpnn_use_side_chain_context 1 --bias_AA "Y:2.0,M:-5.0,C:-5.0" --pack_side_chains 1 --number_of_packs_per_design 4 --pack_with_ligand_context 1 --fixed_residues "{residues}" --repack_everything 0 --batch_size 4 --number_of_batches 5 --checkpoint_path_sc "./model_params/ligandmpnn_sc_v_32_002_16.pt"
"""
    return script.strip() 


def main():
    if len(sys.argv) < 3:
        print("Usage: python ICtoLMPNN.py <protein_PDB_L_dir> <output_dir_fastas>")
        #ex. python ICtoLMPNN.py C:/Users/oliwi/combs/Code/UpdatedIC/Pose31_N2C8O2 C:/Users/oliwi/combs/Code/LigandMPNN/Pose31_N2C802
        #windows user: replace \modified with /modified
        #
        # 
        ##NOTE CHANGE SCRIPT PROCESSING TO BE THE CORRECT OUTPUT DIRECTORY
        return  # Exit the function if there are not enough arguments
    
    protein_L_dir = sys.argv[1]
    if not os.path.isdir(protein_L_dir):
        print(f"Input directory for protein PDBs '{protein_L_dir}' is not a directory.")
        return  # Exit the function if the directory is not valid
    print("protein_pdb_dir is a valid directory")

    fasta_dir = sys.argv[2]
    os.makedirs(fasta_dir, exist_ok=True)
    print(fasta_dir)
    files1 = os.listdir(protein_L_dir)
    counter = 0
    for file in files1:
        counter = counter+1
        if not file.endswith('.pdb'):
            print(f"Skipping non-PDB file: {os.path.basename(file)}")
            continue
        protein_path = os.path.join(protein_L_dir, os.path.basename(file))
        fasta_output_dir = fasta_dir
        os.makedirs(fasta_output_dir, exist_ok=True)
        residues, num_residues = selectE(protein_path, '38E')
        # if num_residues < 1:
        #     continue
        script = ScriptProcessing(protein_path, fasta_output_dir,residues)
        protein_script = os.path.join(fasta_output_dir, f"script{0}.sh")
        print(protein_script)
        with open(protein_script, 'a') as log:
            log.write(script+ '\n')
        print(f"Processing for protein: {os.path.basename(protein_path)}")
        # try:
        #     subprocess.run(protein_script, shell=True, check=True)
        # except subprocess.CalledProcessError as e:
        #     print(f"Error occurred while executing script: {e}")

if __name__ == "__main__":
    main()

# python3 C:/Users/oliwi/LigandMPNN/run.py \
#         --model_type "ligand_mpnn" \
#         --seed 111 \
#         --pdb_path "./inputs/modified_structure_resi_39_GLY_to_ASN_run_0.pdb" \
#         --out_folder "./outputs/sc_fixed_residues" \
#         --pack_side_chains 1 \
#         --number_of_packs_per_design 4 \
#         --pack_with_ligand_context 1 \
#         --fixed_residues "A39 A108" \
#         --repack_everything 0
#         --batch_size 2
#         --number_of_batches 5