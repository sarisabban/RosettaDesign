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

`sudo apt install DSSP`

2. It goes without saying that you need to download and compile [PyRosetta](http://www.pyrosetta.org/) to use this script.

3. You will also need some python libraries that are used within the script, install them using the following commands:

    * If you do not have pip that installs python libraries, install it using this command:
    
    `sudo apt install python3-pip`

    * Install the python libraries
    
    `sudo python3 -m pip install numpy biopython`

## How to use
1. This script and the desired structure to design (FILENAME.pdb) should be in the same working directory.
2. The .pdb file should be cleaned for Rosetta.
3. Run using the following command:

`python3 RosettaDesign.py FILENAME.pdb`

4. The computation takes around 90-120 minutes depending on the protein's size.
5. The script outputs only the final designed .pdb structure.
6. Use [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) to predict the fold of this new design.
