#Monte-Carlo Version/REcursive LMPNNMC.py

#Basically the first input in the fasta file is the basis, then find the best ligand_confidence by looking at the sequence predictability and whichever results in the highest confidence- use that.

import os
import numpy as np
import re
import random
import sys
from pymol import cmd
import matplotlib.pyplot as plt
import subprocess

def selectE(file_path, ligandName):
    print(file_path)
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.select("ligand", f"resname {ligandName}")
    total_residues = cmd.select('temp_selection', 'chain A and name CA')
    print(total_residues)
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
    return residues, num_residues, total_residues


def monte_carlo_ligand_mpnn(sequence, masked):
    masked.sort()
    print("Initial priority queue:", masked)
    iterations = []
    conserved_elements = []
    i = 0
    while True:
        if i >= len(masked) - 1:  # Termination condition based on sequence length
            break
        # Check condition: subsequent number less than the previous
        if masked[i+1] < masked[i]:
            break
        # Simulate sequence conservation: Use current sub-sequence
        conserved_elements.append(masked[i])
        iterations.append(i)

        print(f"Iteration {i}: {masked[i]} conserved.")
        i += 1
        
def parse_data(file_path, masklist): 
    ligand_confidence = []
    fastas = [] #Sequences
    sequence_conservation = [] #Where the sequence is conserved- used to find next masked residues 
    # ligand_confidence.append(ligand_conf_0) #the first ligand is just the base sequence with what has already been masked 
    pattern = r"ligand_confidence=([\d]+\.[\d]{1,4})" #finding the ligand confidence number
    with open(file_path, "r") as file:
        lines = file.readlines()
        first_line = list(lines[1].strip()) #base sequence - full of Gs or As
        counter = 0
        print(first_line)
        for resn in first_line:
            if resn != 'G' and resn != 'A':
                masklist[counter] = resn
            counter = counter + 1
        # fastas.append(first_line) #base sequence
        
        for i in range(2,len(lines),2): #goes through all the odd lines and appends the ligand confidence.
            match = re.search(pattern, lines[i])
            ligand_confidence.append(float(match[1].strip()))
        for i in range(3,len(lines),2): #goes through all the even lines and appends the sequences.
            fastas.append(lines[i].strip())
        
        for i in range(len(fastas[0])): #positions in the sequence, for each residue
            if masklist[i] is not None:
                sequence_conservation.append([{masklist[i] : sum(float(value) for value in ligand_confidence)}])
                continue
            column_dict = {}
            for j in range(len(fastas)): #sequence num
                amino_acid = (fastas[j])[i]
                confidence = ligand_confidence[j]
                #adding aa confidence to dictionary
                if amino_acid in column_dict:
                    column_dict[amino_acid] +=confidence
                else:
                    column_dict[amino_acid] = confidence
                
            # Convert the column dictionary into a list of dictionaries (for each amino acid)
            column_list = [{key: value} for key, value in column_dict.items()]

            # Sort the column_list based on the value in descending order
            column_list.sort(key=lambda x: list(x.values())[0], reverse=True)
            print(column_list)
            sequence_conservation.append(column_list)

        total_ligand_confidence = sum(float(value) for value in ligand_confidence)
        print("TOTAL LIGAND CONFIDENCE")
        print(total_ligand_confidence)
        print(f"SEQUENCE CONSERVATION: {sequence_conservation}")
        counter = 0
        print(len(sequence_conservation))
        for row in sequence_conservation: #these rows should be in order
            print(counter)
            if masklist[counter] is not None:
                print("IS MASKED")
                counter = counter + 1
                continue
            # Extract the value of the first amino acid in each row (if the row is not empty)
            # Check rows to see why some are lists and some are dictionaries
            first_amino_acid_value = 0
            if isinstance(row, list) and isinstance(row[0], dict):
                first_amino_acid_key, first_amino_acid_value =  next(iter(row[0].items()))
            print(f"AMINO ACID CONFIDENCE: {first_amino_acid_value}")
            ran = random.random()
            print(f"RANDOM: {ran}")
            if first_amino_acid_value > (0.9 * total_ligand_confidence) and ran > 0.5:
                print("GREATER CONFIDENCE: ")
                print(f"{list(row[0].keys())[0]}: {first_amino_acid_value}" )
                print(f"{0.8*total_ligand_confidence} : {0.8*total_ligand_confidence < first_amino_acid_value}")
                masklist[counter] = list(row[0].keys())[0]
            counter +=1
        return masklist
                
        
# def output_data(file_path, ligand_conf_0, masklist): 
#     for script in os.listdir(seqs_dir):
#         script_path = os.path.join(outputdir, seqs_dir, script)
#         with open(script_path, 'r') as f_in, open(output, 'a') as f_out:
#             content = ""
#             for _ in range(2):
#                 next(f_in)
#             counter = 0
#             for line in f_in:
#                 if counter % 2 == 0:
#                     name = line.split()
#                     f_out.write((name[0])[:-1]+f"_{int(counter/2)}\n")
#                 else:
#                     f_out.write(line)
#                 counter = counter+1
#             f_out.write("\n")
      
def windows_to_wsl_path(windows_path):
    # Ensure the path uses forward slashes
    wsl_path = windows_path.replace("\\", "/")
    
    # Convert drive letter (e.g., C:) to /mnt/c/
    if wsl_path[1:3] == ":/":
        wsl_path = "/mnt/" + wsl_path[0].lower() + wsl_path[2:]
    
    return wsl_path
      
def ScriptProcessing(pdb_path, output_dir, biasMasking_filepath, residues, counter):
    print(f"SCRIPT PROCESSING: {counter}")
    wsl_pdb_path = windows_to_wsl_path(pdb_path)
    output_path = os.path.join(output_dir, os.path.splitext(os.path.basename(pdb_path))[0], f'{counter}').replace('\\\\', '/')
    output_path = output_path.replace('\\', '/')
    wsl_output_path = windows_to_wsl_path(output_path)
    biasMasking_filepath.replace('\\\\', '/')
    biasMasking_filepath.replace('\\', '/')
    wsl_biasMasking_filepath = windows_to_wsl_path(biasMasking_filepath)    
    if counter != 0:
        script = f"""
#!/bin/sh
python3 /mnt/c/Users/oliwi/LigandMPNN/run.py --model_type "ligand_mpnn" --pdb_path "{wsl_pdb_path}" --out_folder "{wsl_output_path}" --bias_AA_per_residue "{wsl_biasMasking_filepath}" --ligand_mpnn_use_side_chain_context 1 --bias_AA "M:-5.0,C:-5.0" --pack_side_chains 1 --number_of_packs_per_design 2 --pack_with_ligand_context 1 --fixed_residues "{residues}" --repack_everything 0 --batch_size 10 --number_of_batches 10 --checkpoint_path_sc "/mnt/c/Users/oliwi/LigandMPNN/model_params/ligandmpnn_sc_v_32_002_16.pt"
"""
    else:
        script = f"""
#!/bin/sh
python3 /mnt/c/Users/oliwi/LigandMPNN/run.py --model_type "ligand_mpnn" --pdb_path "{wsl_pdb_path}" --out_folder "{wsl_output_path}" --ligand_mpnn_use_side_chain_context 1 --bias_AA "M:-5.0,C:-5.0" --pack_side_chains 1 --number_of_packs_per_design 2 --pack_with_ligand_context 1 --fixed_residues "{residues}" --repack_everything 0 --batch_size 10 --number_of_batches 10 --checkpoint_path_sc "/mnt/c/Users/oliwi/LigandMPNN/model_params/ligandmpnn_sc_v_32_002_16.pt"
"""
    return script.strip() 
            
            
def main():
    if len(sys.argv) < 3:
        print("Usage: python LMPNNMC.py <protein_PDB_dir> <output_dir_name>")
        #ex. python LMPNNMC.py C:/Users/oliwi/combs/Code/UpdatedIC/Pose31_N2C8O2 C:/Users/oliwi/combs/Code/LigandMPNN/Pose31_N2C802
        #windows user: replace \modified with /modified
        #
        # 
        ##NOTE CHANGE SCRIPT PROCESSING TO BE THE CORRECT OUTPUT DIRECTORY
        return  # Exit the function if there are not enough arguments
    
    protein_dir = sys.argv[1]
    if not os.path.isdir(protein_dir):
        print(f"Input directory for protein PDBs '{protein_dir}' is not a directory.")
        return  # Exit the function if the directory is not valid
    print("protein_pdb_dir is a valid directory")
    
    output_dir = os.path.join('C:/Users/oliwi/combs/Code/LigandMPNN', sys.argv[2]).replace('\\\\', '/')
    output_dir = output_dir.replace('\\', '/')
    if os.path.exists(output_dir):
        print(f"Directory {output_dir} already exists. Exiting program. Possibly choose another output directory name")
        sys.exit()  # Exit the program
    os.makedirs(output_dir)
    print(output_dir)
    listPDBS = os.listdir(protein_dir)
    counter = 0
    for pdb in listPDBS:
        counter = counter+1
        if not pdb.endswith('.pdb'):
            print(f"Skipping non-PDB file: {os.path.basename(pdb)}")
            continue
        pdb_path = os.path.join(protein_dir, os.path.basename(pdb)).replace('\\\\', '/')
        pdb_path = pdb_path.replace('\\', '/')
        residues, num_residues, total_residues = selectE(pdb_path, '38E')
        biasMasking_path = os.path.join(output_dir, "biasMasking_0.json") #this is empty for first run, will be filled in subsequent ones.
        open(biasMasking_path, 'w').close()
        script = ScriptProcessing(pdb_path, output_dir, biasMasking_path, residues, 0)
        os.makedirs(os.path.join(output_dir, os.path.splitext(os.path.basename(pdb))[0]))
        os.makedirs(os.path.join(output_dir, os.path.splitext(os.path.basename(pdb))[0], '0')) 
        script_path = os.path.join(output_dir, os.path.splitext(os.path.basename(pdb))[0], '0', f"script_{os.path.splitext(os.path.basename(pdb))[0]}.sh").replace('\\\\', '/')
        script_path = script_path.replace('\\', '/')
        print(script_path)
        with open(script_path, 'a') as log:
            log.write(script+ '\n')
        print(f"Processing for protein: {os.path.splitext(os.path.basename(pdb))[0]}")
        wsl_script_path = windows_to_wsl_path(script_path)
        if not os.path.exists(script_path):
            print("SCRIPT DNE ON WINDOWS")
        else:
            print("SCRIPT EXISTS ON WINDOWS")
        print(wsl_script_path)
        try:
            last_modified_time = os.path.getmtime(script_path)
            # command = ['conda', 'run', '-n', 'ligandmpnn_env', 'wsl', '--exec', '/usr/bin/bash', wsl_script_path]
            # command = f"conda activate ligandmpnn_env && if [ -f {wsl_script_path} ]; then echo 'FILE FOUND'; wsl dos2unix {wsl_script_path} && wsl bash {wsl_script_path}"
            # command = f"conda activate ligandmpnn_env && wsl bash -c if [ -f {wsl_script_path} ]; then dos2unix {wsl_script_path} && cat /mnt/c/Users/oliwi/combs/Code/LigandMPNN/Yang5ExampleMC/3_pose2_no_CG/0/script_3_pose2_no_CG.sh && bash {wsl_script_path}; else echo \"File not found\"; fi'"
            command = (
    f"conda activate ligandmpnn_env && wsl bash -c "
    f"\"cd /mnt/c/Users/oliwi/LigandMPNN && "
    f"if [ -f '{wsl_script_path}' ]; then "
    f"dos2unix '{wsl_script_path}' && "
    f"cat '{wsl_script_path}' && "
    f"bash '{wsl_script_path}'; "
    f"else echo 'File not found'; fi\""
)
            result = subprocess.run(command, capture_output=True, shell=True)
            # subprocess.run(['wsl', 'dos2unix', f'{wsl_script_path}'], check=True)
            # command = [
            #     'wsl', 'bash', '-c',
            #     'source ~/miniconda3/etc/profile.d/conda.sh && conda activate ligandmpnn_env && bash ' + wsl_script_path
            # ]
            #result = subprocess.run(command, capture_output=True)
            print(f"Exit Code: {result.returncode}")
            print(f"STDOUT: {result.stdout.decode()}")
            print(f"STDERR: {result.stderr.decode()}")
            print(f"Successfully executed script using WSL: {script_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while executing script in WSL: {e}")
            sys.exit()
        print("FIRST RUN COMPLETE")
        
        masklist = [None] * total_residues
        for i in range(1,11):
            print(f"On run {i}")
            sequence_path = os.path.join(output_dir, os.path.splitext(os.path.basename(pdb))[0], f'{i-1}', 'seqs', f'{os.path.splitext(os.path.basename(pdb))[0]}.fa').replace('\\\\', '/') #output location for previous sequences.
            sequence_path = sequence_path.replace('\\', '/')
            masklist = parse_data(sequence_path, masklist)
            biasMasking_path = os.path.join(output_dir, f"biasMasking_{i}.json")
            with open(biasMasking_path, 'w') as f:
                f.write("{\n")
                count = 0
                print(masklist)
                for j in range(len(masklist)):
                    if masklist[j] is not None:  # Ensure the value exists
                        print(f"MASKLIST AT: {j+1}: {masklist[j]}")
                        if count == 0:
                            f.write(f'"A{j+1}": {{"{masklist[j]}": 2.0}}')  # Corrected string interpolation
                        else:
                            f.write(f',\n"A{j+1}": {{"{masklist[j]}": 2.0}}')  # Add commas for subsequent items
                        count += 1
                if count > 0:
                    f.write("\n}")  # Close the JSON object if content was written
                else:
                    f.write("}\n")  # Close the JSON if no content is added

            script = ScriptProcessing(pdb_path, output_dir, biasMasking_path, residues, i)
            os.makedirs(os.path.join(output_dir, os.path.splitext(os.path.basename(pdb))[0], f'{i}'))
            script_path = os.path.join(output_dir, os.path.splitext(os.path.basename(pdb))[0], f'{i}', f"script_{os.path.splitext(os.path.basename(pdb))[0]}.sh").replace('\\\\', '/')
            script_path = script_path.replace('\\', '/')
            print(script_path)
            with open(script_path, 'a') as log:
                log.write(script+ '\n')
            print(f"Processing for protein: {os.path.splitext(os.path.basename(pdb))[0]}")
            wsl_script_path = windows_to_wsl_path(script_path)
            try:
                # os.chmod(wsl_script_path, 0o755)  # Grant execute permissions (rwxr-xr-x)
                # command = ['wsl', 'bash', wsl_script_path]
                command = (
    f"conda activate ligandmpnn_env && wsl bash -c "
    f"\"cd /mnt/c/Users/oliwi/LigandMPNN && "
    f"if [ -f '{wsl_script_path}' ]; then "
    f"dos2unix '{wsl_script_path}' && "
    f"cat '{wsl_script_path}' && "
    f"bash '{wsl_script_path}'; "
    f"else echo 'File not found'; fi\""
)
                result = subprocess.run(command, capture_output=True, shell=True)
                print(f"Exit Code: {result.returncode}")
                print(f"STDOUT: {result.stdout.decode()}")
                print(f"STDERR: {result.stderr.decode()}")
                print(f"Successfully executed script using WSL: {script_path}")
                # subprocess.run(command, check=True)
                print(f"Successfully executed script using WSL: {script_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while executing script in WSL: {e}")
            print(f"{i} RUN COMPLETE")
            
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