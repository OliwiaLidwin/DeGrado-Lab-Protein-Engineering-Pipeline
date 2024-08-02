from pymol import cmd

# 
# Used in PlanarDistancy.py.
# Finds all the atoms of interest and gets their coordinates.
# 

def selectE(file_path):
    cmd.reinitialize()
    cmd.load(file_path, "structure")
    cmd.select("ligand", "resname 38E")
    cmd.select("O2_atom", "ligand and name O2")
    cmd.delete("hbonding_oxygens")
    cmd.select("all_oxygens", "elem O")
    
    # print(f"Number of atoms near N: {cmd.count_atoms('all_oxygens')}")
   
    cmd.select("wanted_oxygen", "not resname 38E and not resn SER and (elem O or elem N or elem OH) within 4.4 of O2_atom")
    # cmd.select("wanted_oxygen", "elem O within 3.4 of N_atom and not (resn 38E)")
    NumberofOxygen = cmd.count_atoms('wanted_oxygen')
    
    if NumberofOxygen>1:
        print(f"Number of atoms near N: {NumberofOxygen}")

    cmd.select("C6_atom", "ligand and name C6")
    cmd.select("C7_atom", "ligand and name C7")
    N1_coords = cmd.get_coords("O2_atom")[0].tolist()
    C6_coords = cmd.get_coords("C6_atom")[0].tolist()
    C7_coords = cmd.get_coords("C7_atom")[0].tolist()
    oxygen_atoms = []
    O_coords = []
    exclusion_select=[]
    ran = False

    for run in range(NumberofOxygen):
        ran = True
        if cmd.count_atoms('wanted_oxygen') == 0:
            break        
        O_coords = cmd.get_coords("wanted_oxygen")[0].tolist()
        # print(O_coords)
        oxygen_atoms.append(O_coords)
        # print(oxygen_atoms)
        x,y,z = oxygen_atoms[-1]
        exclusion_select.append(f"(x={x} and y={y} and z={z})")
        exclusion_select_command = " or ".join(exclusion_select)
        cmd.select("wanted_oxygen", f"not ({exclusion_select_command}) and not resname 38E and not resn SER and (elem O or elem N or elem OH) within 3.4 of O2_atom")
        print("appending")
        print(cmd.count_atoms("wanted_oxygen"))
    return N1_coords, C6_coords, C7_coords, oxygen_atoms
