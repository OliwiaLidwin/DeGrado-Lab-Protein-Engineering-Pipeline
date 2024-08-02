import sys
import shutil
import os
import subprocess
from findHBondResi import selectE, numberOfResidues, hasLigand

def scriptWrite(protein_pdbs, protein_pdbs_woLigand, output_dir_scripts, output_dir_proteinPDBS, ligandName):
    apiKey=input("api.register_key: ")
    numOutputs=input("How many outputs? ")
    files1 = os.listdir(protein_pdbs)
    files2 = os.listdir(protein_pdbs_woLigand)

    files1.sort()
    files2.sort()
    
    for file1, file2 in zip(files1, files2):
        if not file1.endswith('.pdb'):
            print(f"Skipping non-PDB file: {os.path.basename(file1)}")
            continue
        file_path1 = os.path.join(protein_pdbs, file1)
        file_path2 = os.path.join(protein_pdbs_woLigand, file2)
        chroma_script_path = os.path.join(output_dir_scripts, os.path.splitext(os.path.basename(file1))[0] + ".py")
        # Make the output directory for the specific protein
        protein_wn_output_direcotry = os.path.join(output_dir_proteinPDBS, os.path.basename(file1))
        os.makedirs(protein_wn_output_direcotry, exist_ok=True)
        # Make a copy of the protein pdb to open it
        pdb_copy = os.path.join(protein_wn_output_direcotry, os.path.basename(file1))
        shutil.copy(file_path1, pdb_copy)
        with open(pdb_copy, 'r') as f:
            pdb_contents = f.read()
            
        # dir_protein_output_path = os.path.join(output_dir_proteinPDBS, os.path.basename(file), f"{counter}")
        if hasLigand(pdb_copy, ligandName):
            print(f"Protein {os.path.basename(file1)} has {ligandName}")
        fixResidues=selectE(pdb_copy, ligandName)
        fixResidues = [int(x) for x in fixResidues]
        proteinLen = numberOfResidues(pdb_copy, ligandName)

        script = f"""import os\nimport torch\nimport chroma\nfrom chroma import Chroma, Protein, api\nfrom chroma.layers.structure.conditioners import SubsequenceConditioner\nfrom chroma.models import graph_classifier, procap, Chroma\napi.register_key('{apiKey}')\nmask=torch.zeros(1,{proteinLen}, device='cuda')\nfixed_residue = {fixResidues}\nmask=torch.ones(1,126, device='cuda')\nmask[0,torch.tensor(fixed_residue)-1]=0\nprotein=Protein.from_PDB('{os.path.basename(file2)}', device='cuda')\nchroma_instance=Chroma()\nfor i in range(10):\n\tprotein1=chroma_instance.design(protein,design_selection=mask)\n\tprotein1.to('{os.path.basename(file2)}')\n\telbo=chroma_instance.score(protein1)['elbo'].score \n\toutput_file=f"{os.path.basename(file2)}_{{i}}_{{elbo}}.pdb"\n\tprotein1.to(os.path.join(output_file))
        
        \nchroma_instance=Chroma()\nnoiser=chroma.layers.structure.diffusion.DiffusionChainCov()\nfor i in range({numOutputs}):\n\tX,C,S=protein.to_XCS()\n\tX_noised=noiser.forward(X, C, 0.5).detach()\n\tnoised_protein = protein.from_XCS(X_noised, C, S)\n\tconditioner = SubsequenceConditioner(chroma_instance.design_network, protein, mask)\n\tprotein1=chroma_instance.sample(protein_init=noised_protein,conditioner=conditioner,tspan=[0.5, 0.001],steps=250,design_selection=mask)\n\toutput_filename = os.path.join('{protein_wn_output_direcotry}', os.path.basename('{file}') + '_' + str(i))\n\telbo = chroma.score(protein1)['elbo'].score\n\tprint(f'sample elbo: {{elbo}}')\n\tprotein1.to(os.path.join(output_filename))"""        
        
        with open(chroma_script_path, 'w') as f:
            f.write(script)
        output = subprocess.run(['python', chroma_script_path], capture_output=True, text=True)
        print(output)
        print(f"Finished chroma run for the protein: {os.path.basename(file1)}")

def main():
    if len(sys.argv) < 6:
        print(len(sys.argv))
        for i in range(len(sys.argv)):
            print((sys.argv[i]))

        print("Usage: python inverseCombsToChromaFixed.py <protein_PBD_dir> <protein_PDB_without_ligand> <output_dir_scripts> <output_dir_protein_PBDs> <ligand_name> ")
        return
    
    protein_pdb_dir = sys.argv[1]
    if not os.path.isdir(protein_pdb_dir):
        print("Input directory for protein PBDs is not a directory.")
        return
    else:
        print("protein_pdb_dir is a valid directory")
        
    protein_pdb_woLigand_dir = sys.argv[2]
    if not os.path.isdir(protein_pdb_woLigand_dir):
        print("Input directory for protein PBDs without ligand is not a directory.")
        return
    else:
        print("protein_pdb_dir without ligand is a valid directory")
    
    output_directory_scripts = sys.argv[3]
    os.makedirs(output_directory_scripts, exist_ok=True)
    output_directory_proteins = sys.argv[4]
    os.makedirs(output_directory_proteins, exist_ok=True)
    scriptWrite(protein_pdb_dir, protein_pdb_woLigand_dir, output_directory_scripts, output_directory_proteins, sys.argv[5])
    print("Finished!")
    
    
if __name__ == "__main__":
    main()
    
    