import math
import sys
import os
import shutil

from selectAtoms import selectE

# 
# Works with selectAtoms.py to find the planar distance created by a ligand based plane to the atom of interest. 
# Creates a new folder in the input directory called all_labeled in which the files are renamed such that the 
# distance for the atom closest to the plane of interest is first in the filename.




def update_and_move_pdb_filename(filename, directory, distances):
    original_path = os.path.join(directory, filename)
    for d in distances:
        filename = f"{os.path.splitext(filename)[0]}_{d}.pdb"
    new_directory = os.path.join(directory, 'all_labeled')
    if not os.path.exists(new_directory):
        os.makedirs(new_directory)
    # print(new_directory)
    shutil.copyfile(original_path, os.path.join(new_directory, filename))
    return filename

def dot_product(v1, v2):
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

def cross_product(v1, v2):
    if len(v1) != 3 or len(v2) != 3:
        raise ValueError("Vectors must be of length 3 to compute cross product.")
    # Compute cross product
    x = v1[1] * v2[2] - v1[2] * v2[1]
    y = v1[2] * v2[0] - v1[0] * v2[2]
    z = v1[0] * v2[1] - v1[1] * v2[0]

    return [x, y, z]

def calculateDistance(co1, co2, co3, co4):  #takes as input the atoms N, C, C, O
    # print(co1)
    # print(co2)
    # print(co3)
    # print(co4)
    d1 = [co2[i] - co1[i] for i in range(3)]
    d2 = [co3[i] - co1[i] for i in range(3)]
    
    #This is fine

    crossprod = cross_product(d1, d2) #this returns the normal vector
    d = (crossprod[0]*-co1[0] + crossprod[1]*-co1[1] + crossprod[2]*-co1[2])
    # print(crossprod)
    # print(d)
    num = abs(dot_product(crossprod, [co4[i] for i in range(3)]) + d) 
    denom = math.sqrt(dot_product(crossprod, crossprod))
    distanceToPlane = num / denom

    return distanceToPlane

def process_files_in_directory(directory):
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        if os.path.isfile(filepath) and filepath.endswith(".pdb"):
            # print(f"Processing file: {filename}")
            coords = selectE(filepath)
            # print("AFTER COORDINATE SELECTION")
            if coords[3]:
                # print(len(coords[3]))
                # print(coords[3])
                distance_values = []
                for i in range(len(coords[3])):
                    # print(len(coords[3][i])//3)
                    distance_values.append(calculateDistance(coords[0], coords[1], coords[2], coords[3][i])) #takes as input the atoms N, C, C, O
                distance_values.sort()
                # print(distance_values)
                new_filename=update_and_move_pdb_filename(filename, directory, distance_values)
                if len(distance_values) > 1:
                    print(f"New file saved as: {new_filename}")
            else:
                continue

def main():
    if len(sys.argv) < 2:
        print("Usage: python PlanarDistance.py <directory_path>")
        return
    directory = sys.argv[1]
    process_files_in_directory(directory)

if __name__ == "__main__":
    main()