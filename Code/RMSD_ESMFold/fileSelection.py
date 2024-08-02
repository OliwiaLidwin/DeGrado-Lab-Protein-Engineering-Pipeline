import os
import re
import sys
import shutil
import matplotlib.pyplot as plt

def get_number_from_filename(filename):
    matches = re.findall(r'(\d+\.\d+)', filename)
    if len(matches) >= 2:
        # Return the last two floating-point numbers as floats
        return float(matches[-2])
    return -1  # Return a large number if the pattern is not found

def create_bar_plot(values):
    values_sorted = sorted(values, key=lambda x: x[0])
    scores = [value for value, filename in values_sorted]
    filenames = [filename for value, filename in values_sorted]
    
    plt.figure(figsize=(12, 8))
    plt.bar(filenames, scores, color='blue')
    plt.xlabel('Filename')
    plt.ylabel('Value')
    plt.title('Values from Filenames')
    plt.xticks(rotation=90)  # Rotate the x labels for better readability
    plt.tight_layout()  # Adjust layout to prevent clipping of tick-labels
    plt.show()
    
def create_bar_plot_20(values):
    values_sorted = sorted(values, key=lambda x: x[0])
    values_sorted = values_sorted[:20]
    scores = [value for value, filename in values_sorted]
    filenames = [filename for value, filename in values_sorted]
    
    plt.figure(figsize=(12, 8))
    plt.bar(filenames, scores, color='blue')
    plt.xlabel('Filename')
    plt.ylabel('Value')
    plt.title('Values from Filenames')
    plt.xticks(rotation=90)  # Rotate the x labels for better readability
    plt.tight_layout()  # Adjust layout to prevent clipping of tick-labels
    plt.show()

def main():
    if len(sys.argv) < 3:
        print("Usage: python fileSelection.py <ESMFoldRMSD_dir> <RenamedOutputDirectory>")
    reference_dir = sys.argv[1]
    if not os.path.isdir(reference_dir):
        print("ESM RMSD(directory) is not a valid directory")
    print(reference_dir)
    renamed_directory = sys.argv[2]    

    values = []
    
    
    listOfESMDirectories = os.listdir(reference_dir)
    for directory in listOfESMDirectories:
        directory_path = os.path.join(reference_dir, directory)
        if os.path.isdir(directory_path):
            listOfFiles = os.listdir(directory_path)
            for file in listOfFiles:
                if file.endswith(".pdb"):
                    value = get_number_from_filename(file)
                    file_path = os.path.join(directory_path, file)
                    output_filename = f"{os.path.basename(directory_path)}_{os.path.basename(file_path)}"
                    values.append((value, output_filename)) 
                    output_path = os.path.join(renamed_directory, output_filename)
                    shutil.copy2(file_path, output_path)
    # listOfESMDirectories = os.listdir(reference_dir)
    # for directory in listOfESMDirectories:
    #     directory_path = os.path.join(reference_dir, directory)
    #     listOfFiles = os.listdir(directory_path)
    #     for file in listOfFiles:
    #         print("FILE PRINTING:")
    #         print(file)
    #         if file.endswith(".pdb"):
    #             value = get_number_from_filename(file)
    #             file_path = os.path.join(directory_path, file)
    #             output_filename = f"{os.path.basename(directory_path)}_{os.path.basename(file_path)[0]}.pdb"
    #             values.append((value, output_filename)) 
    #             os.path.join(renamed_directory, output_filename)
    #             shutil.copy2(file_path, output_filename)
    
    create_bar_plot(values)
    create_bar_plot_20(values)

if __name__ == "__main__":
    main()