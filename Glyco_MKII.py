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
import warnings
warnings.filterwarnings('ignore')

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



def treat_fullfill_dict(protein, THR, carbs, out_infos_buffer, dict_carbs, count):
    """
    Process a single carbohydrate to compute contacts with the protein and update the dictionary.
    """
    carb_positions = np.array([atom.position for atom in carbs.atoms])
    protein_positions = np.array([atom.position for atom in protein.atoms])

    # Compute all pairwise distances between protein and carbohydrate atoms.
    distances = distance_array(protein_positions, carb_positions)

    # Iterate through protein atoms and check for close carbohydrate atoms.
    for i, atom_prot in enumerate(protein.atoms):
        close_carb_indices = np.where(distances[i] < THR)[0]
        if close_carb_indices.size > 0:
            carb_atom = carbs.atoms[close_carb_indices[0]]  # First close carbohydrate atom.
            
            # Write contact info to buffer.
            prot_id = f"{atom_prot.residue.resname}_{atom_prot.residue.resnum}_{atom_prot.segid}"
            carb_info = (f"{prot_id},{carb_atom.segid},{carb_atom.resname},"
                         f"{carb_atom.resid},{carb_atom.type},{count + 1}\n")
            out_infos_buffer.append(carb_info)

            # Update dictionary.
            dict_carbs[carbs.segids[0]].setdefault(prot_id, 0)
            dict_carbs[carbs.segids[0]][prot_id] += 1

def fullfill_dict(THR: float, dict_carbs: dict, SKIP: int):
    """
    Count contacts between carbohydrates and the protein and fulfill the dictionary.

    Parameters:
    -----------
    THR : float
        Threshold distance to consider a contact.
    dict_carbs : dict
        Dictionary to store contact counts.
    SKIP : int
        Number of frames to skip during processing.

    Returns:
    --------
    dict
        Updated dictionary with contact counts.
    """
    # Select carbohydrate atoms (excluding hydrogens).
    input_carbs_list = [u.select_atoms(f"segid {carbs} and not type H") for carbs in dict_carbs.keys()]
    
    # Select C-alpha atoms from the protein.
    protein = u.select_atoms("name CA")

    # Initialize counter for frames and output buffer.
    count = 0
    out_infos_buffer = []

    # Open CSV output file.
    with open("infos_carbos_residue.csv", "w") as out_infos:
        out_infos.write("residue,segid,carbohydrate,carbohydrate_number,group,frame\n")

        # Iterate through trajectory frames.
        for ts in tqdm(u.trajectory):
            if SKIP != 0 and count % SKIP != 0:
                count += 1
                continue

            count += 1
            threads = []
            for carbs in input_carbs_list:
                # Use threading for each carbohydrate set.
                thread = threading.Thread(target=treat_fullfill_dict, args=(protein, THR, carbs, out_infos_buffer, dict_carbs, count))
                threads.append(thread)
                thread.start()

            # Wait for all threads to finish.
            for thread in threads:
                thread.join()

            # Write accumulated results to the file.
            out_infos.writelines(out_infos_buffer)
            out_infos_buffer.clear()

    # Save the contact count dictionary to a CSV file.
    df = pd.DataFrame.from_dict(dict_carbs, orient="index")
    df.to_csv("out_count_carbohydrates.csv")

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

def compute_global_interaction_frequency(interaction_dict, sim_length):
    """
    Compute the global interaction frequency (normalized) for each residue.
    
    Chaque site de glycosylation (par exemple Asn) doit avoir la même valeur de contact, 
    donc nous calculons la fréquence d'interaction uniquement pour chaque site en fonction 
    de ses propres interactions.
    """
    global_frequency = {}  # Initialisation du dictionnaire pour stocker la fréquence globale

    # Itération sur chaque glycan (segment de sucre) et ses interactions
    for carb_data in interaction_dict.values():
        for residue, count in carb_data.items():
            # Normalisation de la fréquence d'interaction pour chaque site de glycosylation
            # La fréquence est calculée en fonction du nombre d'interactions et de la longueur de la simulation
            global_frequency[residue] = global_frequency.get(residue, 0) + count

    # Normalisation de la fréquence (la fréquence doit être entre -1 et 1)
    for residue in global_frequency:
        # On divise par le nombre total de frames de la simulation pour obtenir une fréquence par frame
        global_frequency[residue] = round((global_frequency[residue] / sim_length)*100, 2)
    
    return global_frequency  # Retourne le dictionnaire des fréquences d'interaction globales



def set_global_b_factors(topology, global_data, output_dir):
    """
    Met à jour les b-factors dans le fichier PDB en utilisant les fréquences d'interaction globales.
    """
    output_pdb = f"{output_dir}/global_interaction_frequencies.pdb"  # Définir le nom du fichier de sortie
    parser = PDBParser(QUIET=True)  # Initialiser le parseur PDB
    structure = parser.get_structure("protein", topology)  # Parser la structure PDB

    # Itération sur chaque modèle, chaîne, et résidu dans la structure
    for model in structure:
        for chain in model:
            for residue in chain:
                # Créer un identifiant unique pour chaque résidu (nom, id, et segid)
                res_id = f"{residue.resname}_{residue.get_id()[1]}_{residue.segid}"
                # Récupérer la fréquence d'interaction pour ce résidu, sinon utiliser -1.00 par défaut
                b_factor = global_data.get(res_id, -1.00)  # Si aucune donnée, b-factor par défaut
                # Appliquer la fréquence d'interaction en tant que b-factor pour chaque atome du résidu
                for atom in residue:
                    atom.set_bfactor(b_factor)

    io = PDBIO()  # Initialiser l'objet PDB pour l'écriture
    io.set_structure(structure)  # Définir la structure à sauvegarder
    io.save(output_pdb)  # Sauvegarder la structure mise à jour dans un fichier PDB

def write_log(OUT : str, TOP : str,TRJ :str, THR : int, SKIP : int, full_dict : dict) :
    """Write and print log file.
        
        Parameters
        ----------
            OUT : str
                Path to output files.
            TOP : str
                Path to topology file.
            TRJ : str
                Path to trajectory file.
            THR : int
                Threshold.
            SKIP : int
                Number of frame to skip.
            full_dict : dict
                Dict of count
        
        Returns
        -------
            str
            Write a file.
    """

    log = f"""
        OUTPUT PATH : {OUT}
        TOPOLOGY FILE PATH : {TOP}
        TRAJECTORY FILE PATH : {TRJ}
        THRESHOLD OF CONTACT : {THR} Angstrom
        NUMBER OF FRAME TO SKIP : {SKIP}
        NUMBER OF CARBOHYDATES : {len(full_dict.keys())}
        {[i for i in full_dict.keys()]}
        """
    with open(f"{OUT}/glyco.log", "w") as f :
        f.write(log)
    print(log)
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
    THR = float(THR)
    full_dict = fullfill_dict(THR, dictionnary,SKIP)

    with open('dict_data_julien.pkl', 'wb') as f :
        pkl.dump(full_dict,f)

    #  Frames number.
    if SKIP != 0 :
        full_time = (u.trajectory.totaltime)+1
        full_time = int(full_time/SKIP) + 1
    else :
        full_time = (u.trajectory.totaltime)+1

     # Compute global interaction frequency
    print("Computing global interaction frequency...")
    global_interaction_frequency = compute_global_interaction_frequency(full_dict, full_time)

    # Update PDB file with global frequencies as b-factors
    print("Updating PDB file with global interaction frequencies...")
    set_global_b_factors(args.top, global_interaction_frequency, args.output)


    #  Creation of new structure with new b_factors for each carbohydrates.
    threads = []
    for carb in full_dict.keys() :
        print(f"Treating {carb}...")
        
        #  Create a thread for each carbohydrate to create new structure.
        t = threading.Thread(target=set_new_b_factor,args=(TOP, full_dict, full_time, carb, OUT))
        t.start()
        threads.append(t)
    
    for t in threads :
        t.join()

    #  Write log file.
    write_log(OUT, TOP, TRJ, THR, SKIP, full_dict)
