"""Used to measure presence of carbohydrate on structure.

    Usage
    =====
        ./glyco_diy.py -top <topology_file> -trj <trajectory_file> -output <file_name> -threshold <threshold_interaction> -skip <number_frame_to_skip>
"""

__author__ = "Ravy LEON FOUN LIN"
__date__ = "03 - 10 - 2024"

import argparse
import sys
import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
from math import log
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from Bio.PDB import PDBParser, PDBIO
import pickle as pkl
import threading

parser = argparse.ArgumentParser(prog = 'Glyco if it was actually smart coded', description = 'Compute frequence of interaction by carbohydrate on protein, and set it in new PDB structure as b-factor values.')
parser.add_argument("-top",help="Path to topology file.")
parser.add_argument("-trj",help="Path to trajectory file.")
parser.add_argument("-output",help="Path to output file.")
parser.add_argument("-threshold",help="Threshold of contact.")
parser.add_argument("-skip", help="Set a number if you want to skip a certain number of frames in the simulation (default : 0). Will accelerate the computation.")


args = parser.parse_args()

if args.top == None or args.trj == None or args.output == None or args.threshold == None :
    parser.print_help()
    sys.exit()
else :
    TOP = args.top
    TRJ = args.trj
    OUT = args.output
    THR = args.threshold
    if args.skip != None :
        SKIP = int(args.skip)
    else :
        SKIP = 0

def Universe(TOP : str,TRJ : str) :
    """Set Universe.

        Parameters
        ----------
            TOP : str
                Path to topology file.
            TRJ : str
                Path to trajectory file(s).
    
        Returns
        -------
            Universe (MDAnalysis object.)
    """

    #  Define Universe.
    u = mda.Universe(TOP,TRJ)

    return u


def create_dictionnary(u) :
    """Create Dictionnary of count.

        Parameters
        ----------
            u : MDAnalysis object.
        
        Returns
        -------
            Dictionnary
            Format : {"Carbohydrate_ID" : {'residue':count,...}}
    """

    #  Create empty dictionnary.
    dict_carbs = {}

    #  Fullfill with carbohydrate segid as key and dict of amino acid as values.
    for segment in u.segments :

        #  Get ID treated.
        ID = segment.segid

        #  Check if it is carb
        if ID[:3] == 'CAR' :
            #  Create key-value
            dict_carbs[ID] = {}

    #  Return dict.
    return dict_carbs


def treat_fullfill_dict(protein,THR,carbs,out_infos,dict_carbs,d,count) :
#  Iterate on the different carbohydrates.
    for atom_prot in protein.atoms :
        #  Iterate on each carbohydrate.
        ID_atom_carb = 0
        while d > THR :                    
            #  Compute distance between both atoms.
            d = distance_array(atom_prot.position,carbs[ID_atom_carb].position)[0][0]
            ID_atom_carb += 1
            if ID_atom_carb == len(carbs.atoms) :
                break
        
        #  Check if atom is part of the list.            
        if ID_atom_carb != len(carbs.atoms) :
            #  Check if it under the set threshold.
            if d < THR :
                #  Write the informations file.
                out_infos.write(f"{atom_prot.residue.resname}_{atom_prot.residue.resid}_{atom_prot.segid},{carbs[ID_atom_carb].segid},{carbs[ID_atom_carb].resname},{carbs[ID_atom_carb].resid},{carbs[ID_atom_carb].type},{count+1}\n") 
                #  If it fit in the threshold add to the count.
                if f"{atom_prot.residue.resname}_{atom_prot.residue.resid}_{atom_prot.segid}" not in dict_carbs[carbs.segids[0]].keys() :
                    dict_carbs[carbs.segids[0]][f"{atom_prot.residue.resname}_{atom_prot.residue.resid}_{atom_prot.segid}"] = 1
                else :
                    dict_carbs[carbs.segids[0]][f"{atom_prot.residue.resname}_{atom_prot.residue.resid}_{atom_prot.segid}"] += 1
        #  Then break the for loop, we do not need to count how many atoms of the carbohydrate is in contact. Just if at least one is in contact.
        d = THR +1
    return dict_carbs

def fullfill_dict(THR : str, dict_carbs : dict, SKIP : int) :
    """Count contact of carbohydrate and fullfil dict.

        Parameters
        -----------
            THR : str
                Threshold to count contact.
            dict_carbs : dict
                Dict of count.
            SKIP : int
                Number of frame to skip. (default value : 0)
            
        Returns
        -------
            Dict
            Fullfill dictionnary.
    """

    #  Create list of input for carbohydrates.
    input_carbs_list = [u.select_atoms(f"segid {carbs} and not type H") for carbs in dict_carbs.keys()]
    
    #  Select C-alpha from protein.
    protein = u.select_atoms("name CA")

    #  Count for skipping.
    count = 0

    #  CSV file with contact infos from the carbohydrate.
    out_infos = open("infos_carbos_residue.csv", "w")
    out_infos.write("residue,segid,carbohydrate,carbohydrate_number,group,frame\n")
    #  Distance compared to the threshold.
    d = THR + 1

    #  Iterate on trajectory.
    for ts in tqdm(u.trajectory) :
        
        #  Will treat every SKIP frames.
        if SKIP != 0 :
            if count % int(SKIP) != 0 :
                count += 1
                continue
            else :
                count += 1
        else :
            count += 1
        #  Iterate on protein atoms.
        #  Create a thread for each carbohydrate.
        print(count)
        for carbs in input_carbs_list :
            t = threading.Thread(target=treat_fullfill_dict, args=(protein,THR,carbs,out_infos,dict_carbs,d,count))
            t.start()
            t.join()
    print(dict_carbs)
    #  Save as dataframe.
    df = pd.DataFrame.from_dict(dict_carbs, orient = 'index')
    df.to_csv("out_count_carbohydrates.csv")
    out_infos.close()

    return dict_carbs

def set_new_b_factor(TOP : str, new_b_factors : dict, length_sim : int, carb : str, OUT : str) :
    """New function to set new b-factors values.

        Parameters
        ----------
            TOP : str
                Path to topology file.
            new_b_factors : dict
                Dict which contains count of interaction with carbohydrates?
            length_sim : int
                Duration of the simulation in frame.
            carb : str
                Carbohydrate segid treating.
            OUT : str
                Path to output.

        Returns
        -------
            Write a new file.
    """

    #  Name output.
    output_pdb = f"{OUT}/contact_with_{carb}.pdb"

    #  Set parser.
    parser = PDBParser(QUIET=True)

    #  Parser structure.
    structure = parser.get_structure("protein",TOP)

    #  Full fill dictionnary
    for model in structure :
        for chain in model :
            chain_id = chain.get_id()
            for residue in chain :
                resid = residue.get_id()[1]
                
                resn = residue.resname
                segid = residue.segid
                name = f"{resn}_{resid}_{segid}"
                if name in new_b_factors[carb].keys() :
                    #  Count percentage of interaction through the simulation.
                    new_bfs = (new_b_factors[carb][name]/length_sim)*100
                    
                    #  Round the value.
                    new_bfs = round((new_bfs),2)
                    
                #  If residue never in contact set an impossible value.
                else :
                    new_bfs = -1.00
                #  Replace b-factor value.
                for atom in residue :
                    atom.set_bfactor(new_bfs)
    
    #  Save new structure.
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

    return 0              


if __name__ == "__main__" :
    #  Load universe object...
    print("Load Universe...")
    u = Universe(TOP,TRJ)
    
    #  Create dictionnary for new b_factors.
    print("Create dictionnary...")
    dictionnary = create_dictionnary(u)

    #  Compute contact of carbohydrates with threshold setted.
    print("Fullfill dictionnary...")
    THR = int(THR)
    full_dict = fullfill_dict(THR, dictionnary,SKIP)

    #  Frames number.
    if SKIP != 0 :
        full_time = (u.trajectory.totaltime)+1
        full_time = int(full_time/SKIP) + 1
    else :
        full_time = (u.trajectory.totaltime)+1

    #  Creation of new structure with new b_factors for each carbohydrates?
    for carb in full_dict.keys() :
        print(f"Treating {carb}...")
        
        #  Create a thread for each carbohydrate to create new structure.
        t = threading.Thread(target=set_new_b_factor,args=(TOP, full_dict, full_time, carb, OUT))
        t.start()
        t.join()
