# GlycoFreq
Contains code for the script GlycoFreq. Used to compute from trajectory and topology file, the frequence of residues interacting with carbohydrates.
It produce a csv file with the residue interacting, and the number of frames while they are interacting with the carbohydrates.
But also a PDB file for which underline the residues interacting with the carbohydrates. With the b-factor corresponding to the frequency.

## Settings.

### Environment.
The wide user will create a conda environment in order to safely use this script.

Using the `GlycoFreq.yml` and the following command line, the user can easily create this environment :

`<conda/micromamba> create -n GlycoFreq -f GlycoFreq.yml`

And then :

`<conda/micromamba> activate GlycoFreq`

## Usage.

### CHARMM-GUI.
First of all, you have to treat your topology file with the tool **PDB Reader & Manipulator**.
This will give you a file named `step1_pdbreader.pdb`.
This step is needed in order to 'mark' each carbohydrates with a SEGID, that is used by the script to retrieve the carbohydrates.
If you made your system with **CHARMM-GUI** before, you should already have this file from the .zip they provided. You can also use the file named `step3_input.pdb`, which contains also the water and ions.

### Environment.

Once the environment is activated, here the command line to use **GlycoFreq** :

`./GlycoFreq.py -top <topology_file> -trj <trajectory_file> -output <output_path> -threshold <threshold_interaction> -skip (OPTIONAL) <how_many_frames_to_skip>`
    
### Arguments.
  - **top** : Indicate the path to the topology file (PDB)

  - **trj** : Indicate the path to the trajectory file (XTC,DCD,TRR...).

  - **output** : Indicate the path to the output files. 

  - **threshold** : Indicate the threshold to count a distance as a contact between a heavy atoms of the carbohydrates and the residues.

  - **skip** : Indicate how many frames to skip. **OPTIONAL** 0 by default.


**COOL TRICKS ABOUT TRAJECTORY AND TOPOLOGY FILES** :

  - About the *topology*: To make the computation faster, it is possible to use the topology file with only the protein and carbohydrates. Without the solvant.
    
  - About the *trajectory* : If the topology have no solvant, be sure that you also removed it in the trajectory file (with `gmx trjconv` for example). Also, if the user have replicas, it is possible to concatenate the trajectory files.

## Outputs.

The script produce one *csv* file. 

It contains precise informations about the carbohydrate in interaction with the the residue. Columns of this file => **residue**,**segid**,**carbohydrate**,**carbohydrate_number**,**group**,**frame**.

It produced PDB file of the input protein, but with the b-factor replaced by the percentage of time that a residue was a contact with the query carbohydrate.
One PDB file is produce by carbohydrate. The file is called **contact_with_<segid_of_carbohydrate>.pdb**.

And finally a PDB file of the protein called **global_interaction_frequencies.pdb** is made. This PDB file have the global interactions made by the whole carbohydrates, in the B-factor columns.

And a **log** file.

## Color with Pymol

Once loaded, to show color the protein depending on the frequency of interaction, the following command line can be use in PyMol:

`spectrum b, blue_white_red, selection=polymer`


## Example.

You can try the using the files in the **examples** folder as follow :

`./GlycoFreq.py -top examples/topology_file.pdb -trj examples/trajectory_file.xtc -output . -threshold 8`


And will produce these outputs :

- **contact_with_CARA.pdb** : The PDB file with the proportion of interaction as b-factor.

- **out_count_carbohydrates.csv** : The CSV file which contains the informations of how many time a residue is in contact with the carbohydrate.

- **infos_carbos_residue.csv** : The CSV file which contains the informations about the carbohydrates.

- **glyco.log** : The log file.

These outputs are also in the `examples/results/`.

The PDB file colored by PyMol. It was colored, using the command line shown earlier.
Blue demonstrate a no interactions at all, and going through white to red, demonstrate greater count of contact.
The grey molecule here is the carbohydrate.

<p align="center">
<img src="img/example_structure.gif">
</p>

And the CSV.

**Infos file.**

<p align="center">
<img src="img/example_csv_infos.png">
</p>
