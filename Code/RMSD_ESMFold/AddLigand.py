import os
import re
import sys
import glob
from pymol import cmd

# 
# Adds ligand to ESMFold Output from InverseCombs or LigandMPNN output to be used for MD simulation input.
# Modify lines 14, 23 and 83 accordingly.
# 

def findMatchingFileNames(directory, num, run_number):
    matching_files = []
    pattern = os.path.join(directory, f"m_{num}_*_run_{run_number}_*.pdb")
    print(f"Pattern: {pattern}")

    matching_files.extend(glob.glob(pattern))
    print(matching_files)
    return matching_files

def extractNumberRun(fileName):
    print(fileName)
    pattern = re.compile(r"m_(\d+)_.*_run_(\d+)")
    numbers = (-1,-1)
    match = pattern.search(fileName)
    if match:
        numbers = (int(match.group(1)), int(match.group(2)))
        print(numbers)
    return numbers

def selectE(file_path, reference_path, ligandName, output_dir):
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.load(reference_path, "reference")
    cmd.select("ligand", f"resname {ligandName}")
    cmd.cealign("reference", "structure")
    cmd.sort()
    cmd.rebuild()
    cmd.select("desired_result", "structure or ligand")
    print(cmd.count_atoms("desired_result"))
    cmd.extract("result", "desired_result") #Create new object.
    cmd.remove("reference and not ligand")
    cmd.sort()
    cmd.rebuild()
    print(cmd.count_atoms("result"))
    output_path = os.path.join(output_dir, os.path.basename(file_path))
    cmd.save(output_path, "result")
    print(output_path)
    return output_path
    

def main():
    if len(sys.argv) < 4:
        print("Usage: python AddLigand.py <ESMFOLD_DIRECTORY> <INVERSE_COMBS DIR/ WHERE LIGAND IN PROTEIN> <outputdir>")
        return  # Exit the function if there are not enough arguments
    
    esmfold = sys.argv[1]
    if not os.path.isdir(esmfold):
        print(f"Input directory for protein PDBs '{esmfold}' is not a directory.")
        return  # Exit the function if the directory is not valid
    print("esmfold_dir is a valid directory")

    referencedir = sys.argv[2]
    if not os.path.isdir(referencedir):
        print(f"Input directory for protein PDBs '{referencedir}' is not a directory.")
        return  # Exit the function if the directory is not valid
    print("referencedir is a valid directory")

    outputDir = sys.argv[3]
    os.makedirs(outputDir, exist_ok=True)
    print(outputDir)
    files1 = os.listdir(referencedir)
    for file in files1:
        if not file.endswith('.pdb'):
            print(f"Skipping non-PDB file: {os.path.basename(file)}")
            continue
        protein_path = os.path.join(referencedir, os.path.basename(file))
        runNumber = extractNumberRun(protein_path)
        print(f"{runNumber[0]} run_{runNumber[1]}_")
        MatchingESM = findMatchingFileNames(esmfold, runNumber[0], runNumber[1])
        for j in range(len(MatchingESM)):
            file_path = os.path.join(esmfold, MatchingESM[j])
            selectE(file_path, protein_path, '38E', outputDir)
            
if __name__ == "__main__":
    main()
