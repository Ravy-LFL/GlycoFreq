#!/usr/bin/env python3
"""
Used to measure the presence of carbohydrates on a structure.

Usage:
=====
    ./GlycoFreq.py -top <topology_file> -trj <trajectory_file> -output <output_directory> -threshold <interaction_threshold> -skip <number_of_frames_to_skip>
"""

__author__ = "Ravy LEON FOUN LIN (optimized)"
__date__ = "22-11-2024 (optimized 03-02-2025)"

import argparse
import sys
import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
from Bio.PDB import PDBParser, PDBIO
import threading
import warnings
warnings.filterwarnings('ignore')

# Import the Cython function
import pyximport
pyximport.install(setup_args={"include_dirs": np.get_include()}, language_level=3)
import glyco_distance

# --- Argument parsing ---
parser = argparse.ArgumentParser(prog='GlycoFreq',
                                 description='Compute the frequency of carbohydrate-protein interactions and set them as B-factor values in a new PDB structure.')
parser.add_argument("-top", help="Path to topology file.", required=True)
parser.add_argument("-trj", help="Path to trajectory file.", required=True)
parser.add_argument("-output", help="Output directory.", required=True)
parser.add_argument("-threshold", help="Threshold for contact detection.", required=True)
parser.add_argument("-skip", help="Number of frames to skip (default: 0).", default="0")

args = parser.parse_args()

TOP = args.top
TRJ = args.trj
OUT = args.output
THR = float(args.threshold)
SKIP = int(args.skip)

# --- Main functions ---

def Universe(top: str, trj: str) -> mda.Universe:
    """Creates the Universe object from the topology and trajectory files."""
    return mda.Universe(top, trj)

def create_dictionary(u: mda.Universe) -> dict:
    """Creates an initial dictionary to store contact counts per carbohydrate.
    Format: {"Carbohydrate_ID": { "residue_id": count, ... } }
    """
    return {segment.segid: {} for segment in u.segments if segment.segid.startswith('CAR')}

def treat_fullfill_dict(protein: mda.AtomGroup, THR: float, carbs: mda.AtomGroup, out_infos_buffer: list, dict_carbs: dict, frame_number: int):
    """
    For a given carbohydrate segment, computes contacts with the protein
    using the optimized Cython function and updates the dictionary.
    """
    # Retrieve positions as contiguous arrays
    carb_positions = np.ascontiguousarray(np.array(carbs.positions, dtype=np.float64))
    protein_positions = np.ascontiguousarray(np.array(protein.positions, dtype=np.float64))
    
    # Call the Cython function to get the index of the first contacting carbohydrate atom per protein atom (-1 if no contact)
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
    Iterates over the trajectory and accumulates contact counts per residue for each carbohydrate.
    Only counts the effectively processed frames (i.e., those not skipped).
    """
    # Select carbohydrates (excluding hydrogens)
    input_carbs_list = [u.select_atoms(f"segid {carb} and not type H") for carb in dict_carbs.keys()]
    # Select protein CA atoms
    protein = u.select_atoms("name CA")
    
    frame_number = 0      # Frame counter in the trajectory
    processed_frames = 0  # Number of effectively processed frames
    out_infos_buffer = []
    
    with open("infos_carbos_residue.csv", "w") as out_infos:
        out_infos.write("residue,segid,carbohydrate,carbohydrate_number,group,frame\n")
        for ts in tqdm(u.trajectory, desc="Processing trajectory"):
            frame_number += 1
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
    Computes the global (normalized) interaction frequency for each residue.
    Frequency is defined as (interaction count / number of processed frames) * 100.
    """
    global_frequency = {}
    for carb_data in interaction_dict.values():
        for residue, count in carb_data.items():
            global_frequency[residue] = global_frequency.get(residue, 0) + count
    for residue in global_frequency:
        global_frequency[residue] = round((global_frequency[residue] / processed_frames) * 100, 2)
    return global_frequency

