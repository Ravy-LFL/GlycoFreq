#!/usr/bin/env python3
"""
Used to measure the presence of carbohydrate on a structure.

Usage
=====
    ./GlycoFreq.py -top <topology_file> -trj <trajectory_file> -output <output_directory> -threshold <threshold_interaction> -skip <number_frame_to_skip> -config <path_config_file>
"""

__author__ = "Ravy LEON FOUN LIN"
__date__ = "03-02-2025"
__version__ = "2.0"

import argparse
import sys
import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
from Bio.PDB import PDBParser, PDBIO
import threading
import warnings
import os
warnings.filterwarnings('ignore')

# Import the Cython function
import pyximport
pyximport.install(setup_args={"include_dirs": np.get_include()}, language_level=3)
import glyco_distance

# --- Argument Parsing ---
parser = argparse.ArgumentParser(prog='GlycoFreq',
                                 description='Compute frequency of carbohydrate-protein interactions and set them as B-factor values in a new PDB structure.')
parser.add_argument("-top", help="Path to topology file.", required=True)
parser.add_argument("-trj", help="Path to trajectory file.", required=True)
parser.add_argument("-output", help="Output directory.", required=True)
parser.add_argument("-threshold", help="Threshold of contact.", required=True, default="8.0")
parser.add_argument("-skip", help="Number of frames to skip (default: 0).", default="0")
parser.add_argument("-config",help="Input config file that can be used if no SEGID.", default=False)

args = parser.parse_args()

TOP = args.top
TRJ = args.trj
OUT = args.output
THR = float(args.threshold)
SKIP = int(args.skip)
CONFIG = args.config
# --- Main Functions ---

def Universe(top: str, trj: str) -> mda.Universe:
    """Creates the Universe object from topology and trajectory files."""
    return mda.Universe(top, trj)

def create_dictionary(u: mda.Universe) -> dict:
    """Creates an initial dictionary to store contact counts for each carbohydrate.
    Format: {"Carbohydrate_ID": { "residue_id": count, ... } }
    """
    return {segment.segid: {} for segment in u.segments if segment.segid.startswith('CAR')}

def set_segid(CONFIG: str, TOP: str, TRJ: str) -> mda.Universe:
    """
    Affect segid to the carbohydrates depending on the config file.
    """

    # Load universe.
    u = mda.Universe(TOP, TRJ)

    # Read and treat the configuration file.
    with open(CONFIG, 'r') as f:
        segid_to_residues = {
            line.split(':')[0]: [tuple(res.strip().split('_')) for res in line.split(':')[1].split(',')]
            for line in f
            if line.strip() and not line.startswith('#')
        }

    # Affect the segid.
    for segid, residues in segid_to_residues.items():
        for resn, resi, chain in residues:
            atoms = u.select_atoms(f"resname {resn} and resn {resi} and chainID {chain}")
            # Affect the segid at every atoms selected.
            atoms.segments.segids = [segid] * atoms.n_segments

    return u





def treat_fullfill_dict(protein: mda.AtomGroup, THR: float, carbs: mda.AtomGroup, out_infos_buffer: list, dict_carbs: dict, frame_number: int):
    """
    For a given carbohydrate segment, calculates the contacts with the protein using
    the optimized Cython function and updates the dictionary.
    """
    # Retrieve positions as contiguous arrays
    carb_positions = np.ascontiguousarray(np.array(carbs.positions, dtype=np.float64))
    protein_positions = np.ascontiguousarray(np.array(protein.positions, dtype=np.float64))
    
    # Call the Cython function to obtain, for each protein atom,
    # the index of the first contacting carbohydrate atom (or -1 if no contact).
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
    

def fullfill_dict(THR: float, dict_carbs: dict, SKIP: int, u: mda.Universe, output: str):
    """
    Iterates through the trajectory and accumulates the number of contacts per residue for each carbohydrate.
    Only processes the frames that are not skipped.
    Returns the completed dictionary and the number of processed frames.
    """
    # Select carbohydrates (excluding hydrogens)
    input_carbs_list = [u.select_atoms(f"segid {carb} and not type H") for carb in dict_carbs.keys()]
    # Select CA atoms of the protein
    protein = u.select_atoms("name CA")
    
    frame_number = 0      # Frame counter in the trajectory
    processed_frames = 0  # Number of frames actually processed
    out_infos_buffer = []
    
    with open(f"{output}/infos_carbos_residue.csv", "w") as out_infos:
        out_infos.write("residue,segid,carbohydrate,carbohydrate_number,group,frame\n")
        for ts in tqdm(u.trajectory, desc="Processing trajectory"):
            frame_number += 1
            # Apply skip: process only frames where (frame_number - 1) % SKIP == 0
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
    Calculates the global (normalized) interaction frequency for each residue.
    The frequency is defined as (number of interactions / number of processed frames) * 100.
    """
    global_frequency = {}
    for carb_data in interaction_dict.values():
        for residue, count in carb_data.items():
            global_frequency[residue] = global_frequency.get(residue, 0) + count
    for residue in global_frequency:
        global_frequency[residue] = round((global_frequency[residue] / processed_frames), 2)
    return global_frequency

def write_log(out_dir: str, top: str, trj: str, thr: float, skip: int, full_dict: dict):
    """Writes a log file summarizing the parameters and results."""
    log = f"""
    OUTPUT PATH            : {out_dir}
    TOPOLOGY FILE PATH     : {top}
    TRAJECTORY FILE PATH   : {trj}
    THRESHOLD OF CONTACT   : {thr} Å
    NUMBER OF FRAMES SKIPPED: {skip}
    NUMBER OF CARBOHYDRATES : {len(full_dict.keys())}
    CARBOHYDRATE IDs       : {[i for i in full_dict.keys()]}
    """
    with open(f"{out_dir}/glyco.log", "w") as f:
        f.write(log)
    print(log)
    return 0

def set_new_b_factor(top: str, new_b_factors: dict, sim_length: int, carb: str, out_dir: str):
    """
    Create new PDB File where the B-factor is replaced by the frequency of the carbohydrate with the residue.
    """
    output_pdb = f"{out_dir}/contact_with_{carb}.pdb"

    # Read the input PDB file
    with open(top, 'r') as file:
        lines = file.readlines()

    # Prepare to write the output PDB file
    with open(output_pdb, 'w') as file:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = line[17:20].strip()
                resid = line[22:26].strip()
                segid = line[72:76].strip()  # Assuming segid is in columns 73-76
                name = f"{resname}_{resid}_{segid}"

                # Calculate the new B-factor value
                if name in new_b_factors.get(carb, {}):
                    new_bfs = round(new_b_factors[carb][name] / sim_length, 2)
                else:
                    new_bfs = -1.00

                # Update the B-factor in the line (columns 61-66)
                new_line = f"{line[:60]}{new_bfs:6.2f}{line[66:]}"
                file.write(new_line)
            else:
                file.write(line)


def set_global_b_factors(topology: str, global_data: dict, output_dir: str):
    """
    Update the B-factor with the global coverage of every carbohydrate.
    """
    output_pdb = f"{output_dir}/global_interaction_frequencies.pdb"

    # Read the input PDB file
    with open(topology, 'r') as file:
        lines = file.readlines()

    # Prepare to write the output PDB file
    with open(output_pdb, 'w') as file:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = line[17:20].strip()
                resid = line[22:26].strip()
                segid = line[72:76].strip()  # Assuming segid is in columns 73-76
                res_id = f"{resname}_{resid}_{segid}"

                # Get the new B-factor value
                b_factor = global_data.get(res_id, -1.00)

                # Update the B-factor in the line (columns 61-66)
                new_line = f"{line[:60]}{b_factor:6.2f}{line[66:]}"
                file.write(new_line)
            else:
                file.write(line)


# --- Exécution principale ---
if __name__ == "__main__":

    if not os.path.isdir(OUT):
        os.mkdir(OUT)


    
    if type(CONFIG) == str:
        print("Read config file and Loading Universe...")
        u = set_segid(CONFIG,TOP,TRJ)
    elif CONFIG == False:
        print("Loading Universe...")
        u = Universe(TOP, TRJ)

    print("Creating dictionary...")
    full_dict = create_dictionary(u)

    
    print("Fulfilling dictionary with contact counts...")
    full_dict, processed_frames = fullfill_dict(THR, full_dict, SKIP, u, OUT)
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

