# RosettaDesign
RosettaDesign using PyRosetta

## Decription
This is a python script that allows fixed backbone (fixbb) design of a protein using PyRosetta. Its Design class has three functions:
1. Design.Whole(): Design the whole protein all the onces.
2. Design.Layer(): Designs each protein layer (core, boundry, surface) seperatly.
3. Design.Pack(): Designs each protein layer (core, boundry, surface) seperatly while at the same time designing towards a better packed structure thus resulting in structures with better packed cores. This is the default setting because it is faster and seems to be more accuate (futher testing needed).
4. This script has only been tested in GNU/Linux.

## Requirements
1. You will require DSSP to identify the different protein layers, install using the following command:

`sudo apt install dssp python3-biopython`

2. It goes without saying that you need to download and compile [PyRosetta](http://www.pyrosetta.org/) to use this script.

## How to use
1. This script and the desired structure to design (FILENAME.pdb) should be in the same working directory.
2. The .pdb file should be cleaned for Rosetta.
3. Run using the following command:

`python3 RosettaDesign2.py PROTOCOL FILENAME.pdb`

PROTOCOL will be either fixbb (fix backbone design) or flxbb (flexible backbone design)
FILENAME.pdb the structure in PDB format

4. The computation takes around 1 - 6 hours depending on the protein's size and computer resources.
5. The script outputs only the final designed .pdb structure.
6. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to predict the fold of this new design.
