import sys
import os
import subprocess

from findHBondResi import selectE

def motifscaffolding(protein_pdbs, output_dir):
    files = os.listdir(protein_pdbs)
    protein_pdbs=protein_pdbs+"/"
    dir_makeSec=input("Path of run_inference.py: ")
    model_path=input("Path to Complex_base_ckpt.pt: ")
    for file in files:
        output_prefix = os.path.splitext(file)[0]
        print(output_prefix)

        input_pdb=os.path.join(protein_pdbs, file)
        print("Hydrogen Bonding Results: ")
        print(selectE(input_pdb, '38E'))
        
        protein_design=input("Input contigs: ")
        length = input("Protein length(range): ")
        NumberOfOutputs = input("Number of outputs: ")
        
        script= f"#!/bin/bash \n{dir_makeSec} inference.output_prefix={output_prefix} inference.input_pdb={input_pdb} 'contigmap.contigs=[{protein_design}]' contigmap.length={length} inference.num_designs={NumberOfOutputs} 'inference.ckpt_override_path={model_path}'"
            
        script_path = os.path.join(output_dir, f"{output_prefix}_script.sh")
        with open(script_path, 'w') as f:
            f.write(script)
            
        print(script_path)
        
        os.chmod(script_path, 0o755)
        output = subprocess.run(['bash', script_path], capture_output=True, text=True)
        print(output.stdout)
        if output.stderr:
            print(f"Error for {file}: {output.stderr}")
        
    print(f"Finished making motifscaffolding guideline for {script_path}")
    return 

# def secondaryStructureScript(protein_pdb, output_dir):
#     script_path = os.path.join(output_dir, os.path.basename(protein_pdb))
#     dir_makeSec=input("Path of make_secstruc_adj.py")
    
#     script=f"{dir_makeSec} --pdb_dir={protein_pdb} --out_dir={os.path.dirname(dir_makeSec)}"

#     with open(script_path, 'w') as f:
#         f.write(script)

#     subprocess.run(['python', script_path], capture_output=True, text=True)
    
#     print(f"Finished making secondary structure guideline for {os.path.basename(protein_pdb)}")
#     return os.path.dirname(dir_makeSec)

def main():
    if len(sys.argv) < 3:
        print("Usage: python ICtoRFdiffScriptGen.py <protein_pdb/dir_w_protein_pdbs> <output_dir_for_scripts>")
        return
    
    output_directory = sys.argv[2]
    os.makedirs(output_directory, exist_ok=True)
    
    protein_pdb = sys.argv[1]
    # fileList=[]
    #Make the secondary structure guidelines, put the directory for the output file in a list
    # if os.path.isfile(protein_pdb):
    #     fileList.append(secondaryStructureScript(protein_pdb, output_directory))
    # elif os.path.isdir(protein_pdb):
    #     for file in protein_pdb:
    #         fileList.append(secondaryStructureScript(file, output_directory))
    # else:
    #     return f"{protein_pdb} does not exist or is not a valid path."
    
    if os.path.isdir(protein_pdb):
        motifscaffolding(protein_pdb, output_directory)
    else:
        return f"{protein_pdb} does not exist or is not a valid path."        
    
    
if __name__ == "__main__":
    main()