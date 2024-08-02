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

def get_BB_RMSD(filename):
    matches = re.findall(r'(\d+\.\d+)', filename)
    if len(matches) >= 2:
        # Return the last two floating-point numbers as floats
        return float(matches[-1])
    return -1  # Return a large number if the pattern is not found

def create_bar_plot(values):
    values_sorted = sorted(values, key=lambda x: x[0])
    scores = [value for value, filename in values_sorted]
    filenames = [filename for value, filename in values_sorted]
    
    plt.figure(figsize=(12, 8))
    plt.bar(filenames, scores, color='blue')
    plt.xlabel('Filename')
    plt.ylabel('Side Chain RMSD')
    plt.title('All Pose Side Chain RMSD')
    plt.xticks(rotation=90)  # Rotate the x labels for better readability
    plt.subplots_adjust(bottom=0.6)
    plt.tight_layout()  # Adjust layout to prevent clipping of tick-labels
    plt.savefig(os.path.join('C:/Users/oliwi/combs/Code/RMSD/Best20', 'Side_Chain_RMSD.png'))
    plt.show()
    
def create_bar_plot_20(values, reference_dir, output_dir):
    values_sorted = sorted(values, key=lambda x: x[0])
    values_stored = []
    for i in values_sorted:
        if i[0] < 0.4:
            if get_BB_RMSD(i[1]) < 0.75:
                print(get_BB_RMSD(i[1]))
                values_stored.append(i)
    values_stored = sorted(values_stored, key=lambda x: x[0])
    values_sorted = values_stored[:20]
    scores = [value for value, filename in values_sorted]
    filenames = [filename for value, filename in values_sorted]
    
    for filename in filenames:
        file_path = os.path.join(reference_dir, filename)
        destination = os.path.join(output_dir, filename)
        shutil.copy(file_path, destination)
        
    
    plt.figure(figsize=(12, 8))
    plt.bar(filenames, scores, color='blue')
    plt.xlabel('Filename')
    plt.ylabel('Side Chain RMSD')
    plt.title('20 Best Side Chain RMSD')
    plt.xticks(rotation=90)  # Rotate the x labels for better readability
    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout()  # Adjust layout to prevent clipping of tick-labels
    plt.savefig(os.path.join('C:/Users/oliwi/combs/Code/RMSD/Best20', '20_Best_Side_Chain_RMSD.pdf'))
    plt.show()

def main():
    if len(sys.argv) < 3:
        print("Usage: python RenamedFilesBarPlot.py <RenamedInputDirectory> <OutputDirectory>")
    reference_dir = sys.argv[1]
    if not os.path.isdir(reference_dir):
        print("ESM RMSD(directory) is not a valid directory")
    print(reference_dir)
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)
    values = []
    listOfFiles = os.listdir(reference_dir)
    for file in listOfFiles:
        if file.endswith(".pdb"):
            value = get_number_from_filename(file)
            if value > 0:
                values.append((value, file)) 
    
    create_bar_plot(values)
    create_bar_plot_20(values, reference_dir, output_dir)


if __name__ == "__main__":
    main()