# Glyco-MKII
Contains code for the script Glyco-MKII. Used to compute from trajectory and topology file, the frequence of residues interacting with carbohydrates.
It produce a csv file with the residue interacting, and the number of frames while they are interacting with the carbohydrates.
But also a PDB file for which underline the residues interacting with the carbohydrates. With the b-factor corresponding to the frequency.

## Usage

### Environment.
The wide user will create a conda environment in order to safely use this script.
Using the `Glyco-MKII.yml` and the following command line, the user can easily create this environment :
`<conda/micromamba> create -n Glyco-MKII -f Glyco-MKII.yml`
And then :
`<conda/micromamba> activate Glyco-MKII`

### Inputs.
  - Topology file (PDB,GRO,CIF...).
  - Trajectory file (XTC,DCD,TRR...).
**COOL TRICKS ABOUT THESE FILES** :
  - About the `topology`: To make the computation faster, it is possible to use the topology file with only the protein and carbohydrates. Without the solvant.
  - About the `trajectory` : If the topology have no solvant, be sure that you also removed it in the trajectory file (with `gmx trjconv` for example). Also, if the use have replicas, you concatenate the trajectory files.

