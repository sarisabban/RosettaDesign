# RosettaDesign
RosettaDesign using PyRosetta

## Decription
This is a python script that allows fixed backbone (fixbb) or flexible backbone (flxbb) design of a protein using PyRosetta using any of the following protocols:
1. whole_fixbb: Design the whole protein all the once with a fixed backbone.
2. layer_fixbb: Designs each protein layer (core, boundry, surface) seperatly with a fixed backbone.
3. pack_fixbb: Designs each protein layer (core, boundry, surface) seperatly while at the same time designing towards a better packed structure thus resulting in structures with better packed cores. This is the default setting because it is faster and seems to be more accuate (futher testing needed) with a fixed backbone.
4. pack_flxbb:Designs each protein layer (core, boundry, surface) seperatly while at the same time designing towards a better packed structure thus resulting in structures with better packed cores. This is the default setting because it is faster and seems to be more accuate (futher testing needed) with a flexible backbone.
5. This script has only been tested in GNU/Linux.

## Requirements
1. You will require Biopython and DSSP to identify the different protein layers, install using the following command:

`sudo apt install dssp python3-biopython`

2. You will also need to download and compile [PyRosetta](http://www.pyrosetta.org/).

## How to use
1. It is best if you design proteins that are between 80 - 150 amino acids, any shorter or longer will make it difficult to evaluate the design using the [Rosetta Abinitio Folding Simulation](https://github.com/sarisabban/RosettaAbinitio).
2. This script and the desired structure to design (FILENAME.pdb) should be in the same working directory.
3. The .pdb file should be [cleaned](https://www.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures) for Rosetta/PyRosetta.
4. Run using the following command:

`python3 RosettaDesign.py PROTOCOL FILENAME.pdb`

PROTOCOL will be either fixbb (fix backbone design) or flxbb (flexible backbone design)
FILENAME.pdb the structure in PDB format

5. The computation takes around 1 - 6 hours depending on the protein's size and computer resources.
6. The script outputs only the final designed .pdb structure.
7. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to predict the fold of this new design.
