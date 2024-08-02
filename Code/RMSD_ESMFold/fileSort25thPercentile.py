import os
import re
import sys
import shutil
import numpy as np
import matplotlib.pyplot as plt

def extract_info(filename):
    """Extract relevant parts from the filename and form a new string."""
    base_name = filename.split('.')[0]
    parts = base_name.split('_')
    
    if len(parts) < 7:
        print(f"Warning: Filename {filename} does not have enough parts.")
        return filename

    pose = parts[0]
    resi = parts[1]
    new_residue = parts[4]
    run = '_'.join(parts[5:7])
    
    result = f"{pose}_{resi}_{new_residue}_{run}"
    return result

def get_rmsd_values(filename):
    """Extract SC RMSD and backbone RMSD from filename."""
    matches = re.findall(r'(\d+\.\d+)', filename)
    if len(matches) >= 2:
        sc_rmsd = float(matches[0])  # First float: SC RMSD
        backbone_rmsd = float(matches[1])  # Second float: Backbone RMSD
        return sc_rmsd, backbone_rmsd
    return -1, -1  # Return -1 if pattern is not found or there are not enough numbers

def find_percentile_files(directory, percentile=10):
    """Find files with RMSD values below the specified percentile."""
    all_values = []

    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            sc_rmsd, backbone_rmsd = get_rmsd_values(filename)
            if sc_rmsd != -1 and backbone_rmsd != -1:
                all_values.append((sc_rmsd, backbone_rmsd, filename))

    if not all_values:
        print("No valid files found in the directory.")
        return [], None, None, None, None, None

    sc_rmsd_values = [x[0] for x in all_values]
    backbone_rmsd_values = [x[1] for x in all_values]
    
    if not sc_rmsd_values or not backbone_rmsd_values:
        print("No RMSD values extracted.")
        return [], None, None, None, None, None

    sc_percentile_value = np.percentile(sc_rmsd_values, percentile)
    backbone_percentile_value = np.percentile(backbone_rmsd_values, percentile)

    print(f"{percentile}th percentile for SC RMSD: {sc_percentile_value}")
    print(f"{percentile}th percentile for Backbone RMSD: {backbone_percentile_value}")

    # Identify files below the cutoffs
    below_sc_percentile = [(sc_rmsd, backbone_rmsd, filename) for sc_rmsd, backbone_rmsd, filename in all_values if sc_rmsd < sc_percentile_value]
    below_backbone_percentile = [(sc_rmsd, backbone_rmsd, filename) for sc_rmsd, backbone_rmsd, filename in all_values if backbone_rmsd < backbone_percentile_value]
    both_below = {filename for sc_rmsd, backbone_rmsd, filename in below_sc_percentile} & {filename for sc_rmsd, backbone_rmsd, filename in below_backbone_percentile}
    
    return all_values, below_sc_percentile, below_backbone_percentile, both_below, sc_percentile_value, backbone_percentile_value

def main():
    if len(sys.argv) < 3:
        print("Usage: python fileSort10thPercentile.py <ESMFoldRMSD_dir> <OutputDirForBestResults>")
        return
    
    reference_dir = sys.argv[1]
    if not os.path.isdir(reference_dir):
        print("ESM RMSD directory is not valid.")
        return
    
    output_directory = sys.argv[2]
    os.makedirs(output_directory, exist_ok=True)
    print(f"Output directory created: {output_directory}")
    
    all_files, list_below_sc_percentile, list_below_backbone_percentile, both_below, sc_percentile_value, backbone_percentile_value = find_percentile_files(reference_dir, percentile=10)
    print("BOTH BELOW")
    print(both_below)
    
    if not all_files:
        print("No files found.")
        return
    
    # Renaming and copying files
    renamed_files = []
    for sc_rmsd, backbone_rmsd, filename in list_below_sc_percentile:
        new_filename = extract_info(filename) + ".pdb"
        original_path = os.path.join(reference_dir, filename)
        new_path = os.path.join(output_directory, new_filename)
        shutil.copy2(original_path, new_path)
        renamed_files.append((sc_rmsd, backbone_rmsd, new_filename))
    
    # Debugging output
    print("List of best files after renaming:")
    for sc_rmsd, backbone_rmsd, filename in renamed_files:
        print(f"{filename}: SC RMSD = {sc_rmsd}, Backbone RMSD = {backbone_rmsd}")

    # Plotting SC RMSD vs. Backbone RMSD
    all_sc_rmsd_values = [x[0] for x in all_files]
    all_backbone_rmsd_values = [x[1] for x in all_files]

    plt.figure(figsize=(12, 8))
    plt.scatter(all_sc_rmsd_values, all_backbone_rmsd_values, color='blue', edgecolor='black', label='All Data')

    # Highlight points below Backbone RMSD 10th percentile
    highlighted_sc_rmsd_values = [x[0] for x in list_below_backbone_percentile]
    highlighted_backbone_rmsd_values = [x[1] for x in list_below_backbone_percentile]
    plt.scatter(highlighted_sc_rmsd_values, highlighted_backbone_rmsd_values, color='red', edgecolor='black', label='Below Backbone RMSD 10th Percentile')

    # Highlight points below SC RMSD 10th percentile
    highlighted_sc_rmsd_values = [x[0] for x in list_below_sc_percentile]
    highlighted_backbone_rmsd_values = [x[1] for x in list_below_sc_percentile]
    plt.scatter(highlighted_sc_rmsd_values, highlighted_backbone_rmsd_values, color='green', edgecolor='black', label='Below SC RMSD 10th Percentile')

    # Highlight points below both cutoffs
    both_sc_rmsd_values = [x[0] for x in all_files if x[2] in both_below]
    both_backbone_rmsd_values = [x[1] for x in all_files if x[2] in both_below]
    plt.scatter(both_sc_rmsd_values, both_backbone_rmsd_values, color='purple', edgecolor='black', label='Below Both 10th Percentile')

    # Add percentile lines
    plt.axhline(y=backbone_percentile_value, color='r', linestyle='--', label='Backbone RMSD 10th Percentile')
    plt.axvline(x=sc_percentile_value, color='g', linestyle='--', label='SC RMSD 10th Percentile')
    
    plt.xlabel('SC RMSD')
    plt.ylabel('Backbone RMSD')
    plt.title('SC RMSD vs. Backbone RMSD')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_directory, "Combined_RMSD_Plot.pdf"))
    plt.show()

if __name__ == "__main__":
    main()