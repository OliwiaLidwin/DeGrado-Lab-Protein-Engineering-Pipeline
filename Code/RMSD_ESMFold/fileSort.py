import os
import re
import sys
import heapq
import shutil

def get_number_from_filename(filename):
    matches = re.findall(r'(\d+\.\d+)', filename)
    if len(matches) >= 2:
        # Return the last two floating-point numbers as floats
        return float(matches[-2]), float(matches[-1])
    return -1  # Return a large number if the pattern is not found

def find_lowest_20_files(directory):
    min_heap = []
    
    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            value,value2 = get_number_from_filename(filename)
            
            if value < 0.4:
                if value2 < 0.5:
                # Push tuple (value, filename) onto min_heap
                    heapq.heappush(min_heap, (value, filename))
            
            # # If min_heap exceeds 20 elements, remove the largest element
            # if len(min_heap) > 20:
            #     heapq.heappop(min_heap)

    # Return sorted min_heap which contains the 20 smallest (value, filename) tuples
    return sorted(min_heap)


def main():
    if len(sys.argv) < 3:
        print("Usage: python fileSort.py <ESMFoldRMSD_dir> <OutputDirForBestResults>")
    reference_dir = sys.argv[1]
    if not os.path.isdir(reference_dir):
        print("ESM RMSD(directory) is not a valid directory")
    output_directory = sys.argv[2]
    os.makedirs(output_directory, exist_ok=True)
    print(output_directory)
    listOfBestFiles = find_lowest_20_files(reference_dir)
    
    # Find the 20 files with the lowest values after 'run_'    
    for value, i  in listOfBestFiles:
        print(f"Value: {value}, Filename: {i}")

        original_path = os.path.join(reference_dir, i)
        file_path = os.path.join(output_directory, i)
        shutil.copy2(original_path, file_path)

if __name__ == "__main__":
    main()