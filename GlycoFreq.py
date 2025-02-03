#!/usr/bin/env python3
"""
Used to measure the presence of carbohydrate on a structure.

Usage
=====
    ./optimized_glycofreq.py -top <topology_file> -trj <trajectory_file> -output <output_directory> -threshold <threshold_interaction> -skip <number_frame_to_skip>
"""

__author__ = "Ravy LEON FOUN LIN (optimisé)"
__date__ = "22-11-2024 (optimisé le 03-02-2025)"

import argparse
import sys
import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
from Bio.PDB import PDBParser, PDBIO
import threading
import warnings
warnings.filterwarnings('ignore')

# Import de la fonction Cython
import pyximport
pyximport.install(setup_args={"include_dirs": np.get_include()}, language_level=3)
import glyco_distance

# --- Parsing des arguments ---
parser = argparse.ArgumentParser(prog='GlycoFreq',
                                 description='Compute frequency of carbohydrate-protein interactions and set them as B-factor values in a new PDB structure.')
parser.add_argument("-top", help="Path to topology file.", required=True)
parser.add_argument("-trj", help="Path to trajectory file.", required=True)
parser.add_argument("-output", help="Output directory.", required=True)
parser.add_argument("-threshold", help="Threshold of contact.", required=True)
parser.add_argument("-skip", help="Number of frames to skip (default: 0).", default="0")

args = parser.parse_args()

TOP = args.top
TRJ = args.trj
OUT = args.output
THR = float(args.threshold)
SKIP = int(args.skip)

# --- Fonctions principales ---

def Universe(top: str, trj: str) -> mda.Universe:
    """Crée l'objet Universe à partir des fichiers topologie et trajectoire."""
    return mda.Universe(top, trj)

def create_dictionary(u: mda.Universe) -> dict:
    """Crée un dictionnaire initial pour stocker les comptes de contacts pour chaque sucre.
    Format : {"Carbohydrate_ID": { "residue_id": count, ... } }
    """
    return {segment.segid: {} for segment in u.segments if segment.segid.startswith('CAR')}

def treat_fullfill_dict(protein: mda.AtomGroup, THR: float, carbs: mda.AtomGroup, out_infos_buffer: list, dict_carbs: dict, frame_number: int):
    """
    Pour un segment sucré donné, calcule les contacts avec la protéine en utilisant
    la fonction Cython optimisée et met à jour le dictionnaire.
    """
    # Récupération des positions sous forme de tableaux contigus
    carb_positions = np.ascontiguousarray(np.array(carbs.positions, dtype=np.float64))
    protein_positions = np.ascontiguousarray(np.array(protein.positions, dtype=np.float64))
    
    # Appel à la fonction Cython pour obtenir, pour chaque atome de la protéine,
    # l’indice du premier atome du sucre en contact (ou -1 s'il n'y a pas de contact).
    contacts = glyco_distance.compute_contacts(protein_positions, carb_positions, THR)
    
    for i, atom_prot in enumerate(protein.atoms):
        idx = contacts[i]
        if idx != -1:
            carb_atom = carbs.atoms[idx]
            prot_id = f"{atom_prot.residue.resname}_{atom_prot.residue.resnum}_{atom_prot.segid}"
            carb_info = (f"{prot_id},{carb_atom.segid},{carb_atom.resname},"
                         f"{carb_atom.resid},{carb_atom.type},{frame_number}\n")
            out_infos_buffer.append(carb_info)
            dict_carbs[carbs.segids[0]].setdefault(prot_id, 0)
            dict_carbs[carbs.segids[0]][prot_id] += 1

def fullfill_dict(THR: float, dict_carbs: dict, SKIP: int, u: mda.Universe):
    """
    Parcourt la trajectoire et cumule le nombre de contacts par résidu pour chaque sucre.
    Ne compte que les frames effectivement traitées (c'est-à-dire celles qui ne sont pas ignorées).
    Retourne le dictionnaire complété et le nombre de frames traitées.
    """
    # Sélection des sucres (exclusion des hydrogènes)
    input_carbs_list = [u.select_atoms(f"segid {carb} and not type H") for carb in dict_carbs.keys()]
    # Sélection des atomes CA de la protéine
    protein = u.select_atoms("name CA")
    
    frame_number = 0      # Compteur du numéro de frame dans la trajectoire
    processed_frames = 0  # Nombre de frames réellement traitées
    out_infos_buffer = []
    
    with open("infos_carbos_residue.csv", "w") as out_infos:
        out_infos.write("residue,segid,carbohydrate,carbohydrate_number,group,frame\n")
        for ts in tqdm(u.trajectory, desc="Processing trajectory"):
            frame_number += 1
            # Appliquer le skip : traiter uniquement les frames où (frame_number - 1) % SKIP == 0
            if SKIP and (frame_number - 1) % SKIP != 0:
                continue

            processed_frames += 1
            threads = []
            for carbs in input_carbs_list:
                t = threading.Thread(target=treat_fullfill_dict,
                                     args=(protein, THR, carbs, out_infos_buffer, dict_carbs, processed_frames))
                threads.append(t)
                t.start()
            for t in threads:
                t.join()
            out_infos.writelines(out_infos_buffer)
            out_infos_buffer.clear()
    
    return dict_carbs, processed_frames

def compute_global_interaction_frequency(interaction_dict: dict, processed_frames: int):
    """
    Calcule la fréquence globale (normalisée) d'interaction pour chaque résidu.
    La fréquence est définie comme (nombre d'interactions / nombre de frames traitées) * 100.
    """
    global_frequency = {}
    for carb_data in interaction_dict.values():
        for residue, count in carb_data.items():
            global_frequency[residue] = global_frequency.get(residue, 0) + count
    for residue in global_frequency:
        global_frequency[residue] = round((global_frequency[residue] / processed_frames) * 100, 2)
    return global_frequency

def set_new_b_factor(top: str, new_b_factors: dict, sim_length: int, carb: str, out_dir: str):
    """
    Crée une nouvelle structure PDB en affectant pour chaque résidu le pourcentage d'interaction en tant que B-factor.
    """
    output_pdb = f"{out_dir}/contact_with_{carb}.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", top)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                resid = residue.get_id()[1]
                resn = residue.resname
                segid = residue.segid
                name = f"{resn}_{resid}_{segid}"
                if name in new_b_factors[carb]:
                    new_bfs = round((new_b_factors[carb][name] / sim_length) * 100, 2)
                else:
                    new_bfs = -1.00
                for atom in residue:
                    atom.set_bfactor(new_bfs)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    return 0

def set_global_b_factors(topology: str, global_data: dict, output_dir: str):
    """
    Met à jour les B-factors dans le fichier PDB en utilisant la fréquence d'interaction globale.
    """
    output_pdb = f"{output_dir}/global_interaction_frequencies.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", topology)
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = f"{residue.resname}_{residue.get_id()[1]}_{residue.segid}"
                b_factor = global_data.get(res_id, -1.00)
                for atom in residue:
                    atom.set_bfactor(b_factor)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def write_log(out_dir: str, top: str, trj: str, thr: float, skip: int, full_dict: dict):
    """Écrit un fichier log résumant les paramètres et résultats."""
    log = f"""
OUTPUT PATH            : {out_dir}
TOPOLOGY FILE PATH     : {top}
TRAJECTORY FILE PATH   : {trj}
THRESHOLD OF CONTACT   : {thr} Å
NUMBER OF FRAMES SKIPPED: {skip}
NUMBER OF CARBOHYDATES : {len(full_dict.keys())}
CARBOHYDATES IDs       : {[i for i in full_dict.keys()]}
"""
    with open(f"{out_dir}/glyco.log", "w") as f:
        f.write(log)
    print(log)
    return 0

# --- Exécution principale ---
if __name__ == "__main__":
    print("Loading Universe...")
    u = Universe(TOP, TRJ)
    
    print("Creating dictionary...")
    full_dict = create_dictionary(u)
    
    print("Fulfilling dictionary with contact counts...")
    full_dict, processed_frames = fullfill_dict(THR, full_dict, SKIP, u)
    print(f"Number of processed frames: {processed_frames}")
    
    print("Computing global interaction frequencies...")
    global_interaction_frequency = compute_global_interaction_frequency(full_dict, processed_frames)
    
    print("Updating global B-factors in PDB file...")
    set_global_b_factors(TOP, global_interaction_frequency, OUT)
    
    print("Creating new structures with per-carbohydrate B-factors...")
    threads = []
    for carb in full_dict.keys():
        print(f"Processing carbohydrate {carb}...")
        t = threading.Thread(target=set_new_b_factor, args=(TOP, full_dict, processed_frames, carb, OUT))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    
    write_log(OUT, TOP, TRJ, THR, SKIP, full_dict)

