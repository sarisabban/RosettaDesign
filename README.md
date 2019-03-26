# RosettaDesign
RosettaDesign using PyRosetta

## Decription
This is a python script that allows fixed backbone (fixbb) or flexible backbone (flxbb) design of a protein in PyRosetta.

For fixbb:
1. Relax structure.
2. Fixbb (designs the whole protein keeping the backbone fixed)
3. Relax structure.

For flxbb (recommended):
1. Relax structure.
2. Flxbb (designs the whole protein while allowing for a flexible backbone).
3. Relax designed structure.

This script has only been tested in GNU/Linux.

## Requirements
1. You will need to download and compile [PyRosetta](http://www.pyrosetta.org/).
2. You will also need to install Biopython and dssp to run the BLAST operation after the RosettaDesign protocol, install using the following command:

`sudo apt install python3-biopython dssp`

## How to use
1. It is best if you design proteins with a single chain that is between 80 - 150 amino acids, any shorter or longer will make it difficult to evaluate the design using the [Rosetta Abinitio Folding Simulation](https://github.com/sarisabban/RosettaAbinitio).
2. This script and the desired structure to design (FILENAME.pdb) should be in the same working directory.
3. The .pdb file should be [cleaned](https://www.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures) for Rosetta/PyRosetta.
4. Run using the following command:

`python3 RosettaDesign.py --PROTOCOL FILENAME.pdb`

* PROTOCOL = Will be either --fixbb (for fix backbone design) or --flxbb (for flexible backbone design)
* FILENAME.pdb = The structure to be designed in PDB format

5. The computation takes around ~100 hours (fixbb) or ~200 hours (flxbb) depending on the protein's size and computer resources.
6. The script will output around 100 structures and a score file *.fasc*, then it will auto choose the lowest energy scoring structure and delete all the others, including the scoring file.
7. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to predict the fold of the new desgined lowest scoring structure.
8. NOTE: Some backbones are not possible to simulate their folding (*Abinitio*) before/after design even if the designing step was successful, this is as a result of the fragments used to to predict the structure's fold not being optimal and not an issue with the design algorithm, so far there is no way arround this (this is a limitation in Rosetta/PyRosetta).
