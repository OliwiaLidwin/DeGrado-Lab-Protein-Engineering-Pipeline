import os
import shutil

def match_and_copy_files(dir1, dir2, output_dir):
    """
    Check if files with the same name exist in both directories and copy them to the output directory.
    
    Args:
        dir1 (str): Path to the first directory.
        dir2 (str): Path to the second directory.
        output_dir (str): Path to the directory where matching files will be copied.
        
    Change lines 47 to 49
    """
    # Get the set of filenames from both directories
    files_dir1 = set(os.listdir(dir1))
    files_dir2 = set(os.listdir(dir2))
    
    print(f"Files in {dir1}: {files_dir1}")
    print(f"Files in {dir2}: {files_dir2}")
    
    # Find common filenames
    common_files = files_dir1.intersection(files_dir2)
    
    print(f"Common files: {common_files}")
    
    if not common_files:
        print("No matching files found in both directories.")
        return
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    for filename in common_files:
        # Define paths for source and destination files
        file_path_dir1 = os.path.join(dir1, filename)
        file_path_dir2 = os.path.join(dir2, filename)
        output_file_path = os.path.join(output_dir, filename)
        
        # Copy files from both directories to the output directory
        shutil.copy2(file_path_dir1, output_file_path)
        shutil.copy2(file_path_dir2, output_file_path)
        
        print(f"Copied {filename} to {output_dir}")

# Example usage:
dir1 = r"C:/Users/oliwi/combs/Code/RMSD/esmfold7_22/renamed/renamedWLigand/best_t10/best_SCRMSD"
dir2 = r"C:/Users/oliwi/combs/Code/RMSD/esmfold7_22/renamed/renamedWLigand/best_t10/best_BBRMSD"
output_dir = r"C:/Users/oliwi/combs/Code/RMSD/esmfold7_22/renamed/renamedWLigand/best_t10/best_match"

match_and_copy_files(dir1, dir2, output_dir)
