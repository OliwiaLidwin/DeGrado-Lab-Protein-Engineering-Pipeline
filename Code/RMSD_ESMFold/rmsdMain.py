import os
import sys
import shutil
import glob
import re
from rmsdCalcInitial import RMSD_Calc

# 
# Main file for calculating the backbone and side chain RMSD and renaming the filename with the numbers. 
# Currently the side chain measurement is commented out and only the backbone RMSD is usable. 
# Update lines 19,33 accordingly.
# 


# def findMatchingFileNames(directory, pose_number, num, run_number):
def findMatchingFileNames(directory, pose_number):
    matching_files = []
    # pattern = os.path.join(directory, f"Pose{pose_number}_{num}_*_run_{run_number}_*.pdb")
    pattern = os.path.join(directory, f"07_30_32_902_esmfold_{pose_number}_A.pdb")
    print(f"Pattern: {pattern}")
    matching_files.extend(glob.glob(pattern))
    print(matching_files)
    return matching_files

def extractNumberRun(fileName):
    print(fileName)
    # pattern = re.compile(r"Pose(\d+)_(\d+)_.*_run_(\d+)")
    # numbers = (-1,-1)
    # match = pattern.search(fileName)
    # if match:
    #     numbers = (int(match.group(1)), int(match.group(2)), int(match.group(3)))
    #     print(numbers)
    pattern = re.compile(r"07_30_32_902_esmfold_(\d+).pdb")
    numbers = (-1)
    match = pattern.search(fileName)
    if match:
        numbers = (int(match.group(1)))
        print(numbers)
    return numbers

def main():
    if len(sys.argv) < 4:
        print("Usage: python rmsdMain.py <ICResultPDB(directory)> <dirOfESMFOLD> <OutputDir>")
    
    #Check that everything is valid
    reference_pdb_path = sys.argv[1]
    if not os.path.isdir(reference_pdb_path):
        print("Reference PDB (directory) is not a valid directory")
    result_directory = sys.argv[2]
    if not os.path.isdir(result_directory):
        print("ESMFOLD directory is not a directory")
    output_directory = sys.argv[3]
    os.makedirs(output_directory, exist_ok=True)
    print(output_directory)
    
    #These should be in the same order
    listICresults = os.listdir(reference_pdb_path)
    listESMPDB = os.listdir(result_directory)
    listICresults = sorted(listICresults)
    listESMPDB = sorted(listESMPDB)
        
    for i in range(len(listICresults)):
        if not listICresults[i].endswith(".pdb"):
            continue
        runNumber = extractNumberRun(listICresults[i])
        reference_pdb = os.path.join(reference_pdb_path, listICresults[i])
        # print(f"Pose_{runNumber[0]} {runNumber[1]} run_{runNumber[2]}_")
        print(f"07_30_32_902_esmfold_{runNumber}")

        # MatchingESM = findMatchingFileNames(result_directory, runNumber[0], runNumber[1], runNumber[2])
        MatchingESM = findMatchingFileNames(result_directory, runNumber)

        for j in range(len(MatchingESM)):
            file_path = os.path.join(result_directory, MatchingESM[j])
            # value_sc, value_bb  = RMSD_Calc(reference_pdb, file_path)
            value_bb  = RMSD_Calc(reference_pdb, file_path)
            # print(value_sc)
            print(value_bb)
            print(MatchingESM[j])
            # outputFileName = f"{os.path.splitext(os.path.basename(MatchingESM[j]))[0]}_{value_sc}_{value_bb}.pdb"
            outputFileName = f"{os.path.splitext(os.path.basename(MatchingESM[j]))[0]}_{value_bb[1]}.pdb"
            print(outputFileName)
            outputFile = os.path.join(output_directory, outputFileName)
            shutil.copy2(file_path, outputFile)
            # print(f"Copdied contents from {file_path} with rmsd_sc {value_sc} and rmsd_bb {value_bb} to {outputFile}")
            print(f"Copdied contents from {file_path} with rmsd_bb {value_bb} to {outputFile}")

if __name__ == "__main__":
    main()