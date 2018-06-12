# RosettaDesign
RosettaDesign using PyRosetta

## Decription
This is a python script that allows fixed backbone (fixbb) or flexible backbone (flxbb) design of a protein using PyRosetta using any of the following protocols:
1. fixbb: Designs the whole protein keeping the backbone fixed.
2. flxbb: Designs the whole protein while allowing for a flexible backbone.
3. This script has only been tested in GNU/Linux.

## Requirements
1. You will need to download and compile [PyRosetta](http://www.pyrosetta.org/).
2. You will also need to install Biopython to run the BLAST operation after the RosettaDesign protocol, install using the following command:

`sudo apt install python3-biopython`

## How to use
1. It is best if you design proteins with a single chain that is between 80 - 150 amino acids, any shorter or longer will make it difficult to evaluate the design using the [Rosetta Abinitio Folding Simulation](https://github.com/sarisabban/RosettaAbinitio).
2. This script and the desired structure to design (FILENAME.pdb) should be in the same working directory.
3. The .pdb file should be [cleaned](https://www.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures) for Rosetta/PyRosetta.
4. Run using the following command:

`python3 RosettaDesign.py PROTOCOL FILENAME.pdb`

* PROTOCOL = Will be either fixbb (for fix backbone design) or flxbb (for flexible backbone design)
* FILENAME.pdb = The structure to be designed in PDB format

5. The computation takes around 24-48 hours depending on the protein's size and computer resources.
6. The script will preform 50 relax operations (taking the lowest scoring structure) then 100 design operations (taking the lowest scoring structure).
7. The script outputs only one file which is the final designed structure named *structure.pdb*.
8. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to predict the fold of your new desgined structure.
