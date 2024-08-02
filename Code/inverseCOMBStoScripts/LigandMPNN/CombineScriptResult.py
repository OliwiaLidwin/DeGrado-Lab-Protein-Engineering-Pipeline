import sys
import os

# 
# Groups all LigandMPNN output fasta sequences into one file. 
# Output of this file should be a fasta file to 
# 

def list_directories(directory):
    if os.path.exists(directory):
        directories = [name for name in os.listdir(directory) if os.path.isdir(os.path.join(directory, name))]
        return directories
    else:
        print(f"The directory '{directory}' does not exist.")
        return []


def main():
    if len(sys.argv) < 3:
        print("Usage: python CombineScriptResult.py <../../script0.sh (parent of the parent directory)> <outputdir>")
        return  # Exit the function if there are not enough arguments
    
    dd1scripts = sys.argv[1] #Pose95_N2C8O2
    if not os.path.isdir(dd1scripts):
        print(f"Input ../../script '{dd1scripts}' is not a directory.")
        return  # Exit the function if the directory is not valid
    print(f"{dd1scripts} is a valid directory") #this is good
    
    
    outputdir = sys.argv[2]
    if not os.path.isdir(outputdir):
        print(f"Input '{outputdir}' is not a directory.")
        return  # Exit the function if the directory is not valid
    print(f"{outputdir} is a valid directory") #this is good
    
    output = os.path.join(outputdir, "output.fa")   
    print(output)
    existing_FA = []
    
    for dscripts in os.listdir(dd1scripts): #iteratres through the list of objects
        dscripts_path = os.path.join(outputdir, dd1scripts, dscripts)
        if not os.path.isdir(dscripts_path):
            continue
        print(dscripts)
        
        
        seqs_dir = os.path.join(outputdir, dscripts, "seqs")
        print(seqs_dir)
        os.makedirs(seqs_dir, exist_ok=True)
        
        for script in os.listdir(seqs_dir):
            script_path = os.path.join(outputdir, seqs_dir, script)
            with open(script_path, 'r') as f_in, open(output, 'a') as f_out:
                content = ""
                for _ in range(2):
                    next(f_in)
                counter = 0
                for line in f_in:
                    if counter % 2 == 0:
                        name = line.split()
                        f_out.write((name[0])[:-1]+f"_{int(counter/2)}\n")
                    else:
                        f_out.write(line)
                    counter = counter+1
                f_out.write("\n")
    print("Combination complete")
            
if __name__ == "__main__":
    main()



# >modified_structure_resi_39_GLY_to_GLN_run_3, T=0.1, seed=111, num_res=124, num_ligand_res=28, use_ligand_context=True, ligand_cutoff_distance=8.0, batch_size=2, number_of_batches=5, model_path=./model_params/ligandmpnn_v_32_010_25.pt
# GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGQGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGSGGGGGGGGGGGGGGGGGG
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=1, T=0.1, seed=111, overall_confidence=0.3609, ligand_confidence=0.3725, seq_rec=0.0403
# SLKSRFEELKKVFDKIVETLKKMEKAFEDGDLEEVAKLQAEIKKLVAEEKKLRKKLLEEAKKAGNTEAVKLMEEYEKVREEGVKLGKEADKAFKAGDKAAVAAYLKQSRELLKQCEKLLKEIEKAL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=2, T=0.1, seed=111, overall_confidence=0.3606, ligand_confidence=0.3740, seq_rec=0.0323
# SLKSKFEEMKKVFDKIVATLKEMKKAWEDGDWETVAKLQEEIKKLIAEYKKLLEKLLKEAEKAGLTEAVKLMREYEKARDEAVAIGAKADAAFKAGDEAAVAAYLKQSEELLKKCEKLLKEIEKAL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=3, T=0.1, seed=111, overall_confidence=0.3667, ligand_confidence=0.3868, seq_rec=0.0403
# SLKSLFEEMKKVYDKIVATLEKMEKAFEDGDMEEVEKLQKEIEKLIKEYKKLLEKLLEEAKAAGNEEAVKLMEELKKVRDEGEALGKETIKAIKAGDTAKVKAYLEQSKALRAKCDTLMKEVEKAL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=4, T=0.1, seed=111, overall_confidence=0.3456, ligand_confidence=0.3671, seq_rec=0.0565
# SLKSLFEELLKVYEETVKTLGEMVAAFEAGDMATVAELQARIEKLIAEEKKLLKVLLEQARAAGNTEAVALMEEYEKVRNEGIAIGKAANAALAAGDHAAVAAYIAQSQALLAQCDVLLEEIGKAL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=5, T=0.1, seed=111, overall_confidence=0.3429, ligand_confidence=0.3760, seq_rec=0.0484
# SLESQFEQLQAVYDEIVATLKEMVAAWKAGDRERVAELQDRITGLIAEEKALRAKLLEQARAAGNTEAVALLERYEAVRDEGTAIGAESRAAFEAGDDEAVAAHIKRSEELLAKCQKLLEEIRKAL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=6, T=0.1, seed=111, overall_confidence=0.3570, ligand_confidence=0.3408, seq_rec=0.0403
# SLESLFEEMKKVFDELVATLKEMVAAFKAGDMAKVAELQARIDKLIAELNKLRKKLLALAKAAGNTEAVKLMEEYEKVRNEGVAVGKATDAAIAAGDTAAVEAHIKKSEELLKKCEKLLKEIEEAL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=7, T=0.1, seed=111, overall_confidence=0.3348, ligand_confidence=0.3547, seq_rec=0.0484
# SLKSQFQELKKVFDETVATLGKMEKAFKNGDTATVAKLQKTIDKLIAEEQKLLKKLLEEAKAAGNTEAVTLLEQYEAARNQGVAVGKKSDAAFKAGDTAAVEAYLKQSRTLLAQCQTLLKQIEQSL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=8, T=0.1, seed=111, overall_confidence=0.3372, ligand_confidence=0.3500, seq_rec=0.0484
# SLKSLFEELLKVYEKIVATLGEMKAAFEAGDMDRVRELQEEIVKLIKEAKKLLEKLLAMARAAGNTEAVALMEKYREVREKGEAVGKKSIEAIEAGDWEAVRAYLAQSDELLAQCEKLLEEIEKAL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=9, T=0.1, seed=111, overall_confidence=0.3420, ligand_confidence=0.3579, seq_rec=0.0484
# SLESLFQEMLKVYNEIVKTLKEMEAAFKNGDFAKVKELQAKIEKLIAELKKLEEKLLAEARAAGNTEAVALMEEYRKVREEGEKVGAEAIKAFEAGDWKEYEKLLKKSEELLKKCDTLLKKIGAAL
# >modified_structure_resi_39_GLY_to_GLN_run_3, id=10, T=0.1, seed=111, overall_confidence=0.3650, ligand_confidence=0.3518, seq_rec=0.0484
# SLESKFEELLKVYEKTVETLKKMEEALKNGDMEKVAELQEEIKKLIKEEKKLLEELLEEAREAGLTEAVKLMEEYKKVREKGEKVGEKSVEAFKAGDLAEVEKYLEESKKLLKECEVLLEEIGKAL