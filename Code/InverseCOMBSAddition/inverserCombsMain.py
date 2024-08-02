import sys
import os
import gzip
from inverseCombsToPDB import replaceResidue

# 
# Adds the inverseCombsOutput from .gzip to the specific(COMBS) output file.
# 


def main():
    counter = 0
    if len(sys.argv) < 4:
        print("Usage: python inverserCombsMain.py <protein_pdb> <dir_inverseCombs_output> <output_dir>")
        #ex. python inverserCombsMain.py C:/Users/oliwi/Downloads/All/yang_31_pose101_no_CG.pdb  C:/Users/oliwi/combs/inverseCombs/temp/yang_31_pose101_no_CG_BB_CCO C:/Users/oliwi/combs/Code/outputs/Pose31_N2C8O2
        return
    
    #Check if protein pdb is valid
    protein_pdb = sys.argv[1]
    if not os.path.isfile(protein_pdb):
        print("Protein file is not found")
        sys.exit(1)
    elif not protein_pdb.endswith(".pdb"):
        print("Protein PDB is not a pdb file")
        sys.exit(1)


    dir_inverseCombs_output = sys.argv[2]
    for gzipFilename in os.listdir(dir_inverseCombs_output):
        gzipFilePath = os.path.join(dir_inverseCombs_output, gzipFilename)
        if not gzipFilename.endswith(".gz"):
            print(f"Error: {gzipFilePath} is not a valid gzip file. Skipping.")
            continue
        
        #Extracting pdb from zip file
        with gzip.open(gzipFilePath, 'rb') as gz_file:
            pdb_content = gz_file.read().decode('utf-8')
            
            temp_file_path = os.path.join(dir_inverseCombs_output, gzipFilename[:-3])

            with open(temp_file_path, 'w') as temp_file:
                temp_file.write(pdb_content)

            if not os.path.isfile(temp_file_path):
                print(f"Error: Failed to create temp file {temp_file_path}")
                continue
            
            #Replacing the residue
            output_path = replaceResidue(protein_pdb, temp_file_path, sys.argv[3],counter)
            print(f"Saved the output to: {output_path}")
            counter=counter+1
            
if __name__ == "__main__":
    main()