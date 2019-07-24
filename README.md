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

`sudo apt install python3-biopython dssp gnuplot`

## How to use
1. It is best if you design proteins with a single chain that is between 80 - 150 amino acids, any shorter or longer will make it difficult to evaluate the design using the [Rosetta Abinitio Folding Simulation](https://github.com/sarisabban/RosettaAbinitio).
2. This script and the desired structure to design (FILENAME.pdb) should be in the same working directory.
3. The .pdb file should be [cleaned](https://www.rosettacommons.org/docs/latest/rosetta_basics/preparation/preparing-structures) for Rosetta/PyRosetta.
4. Run using the following command:

`python3 RosettaDesign.py --PROTOCOL FILENAME.pdb`

* PROTOCOL = Will be either --fixbb (for fix backbone design) or --flxbb (for flexible backbone design)
* FILENAME.pdb = The structure to be designed in PDB format

5. The computation takes around ~100 hours (fixbb) or ~300 hours (flxbb) depending on the protein's size and computer resources.
6. The script will output around 100 structures and a score file *.fasc*.
5. Use this [Rosetta Abinitio](https://github.com/sarisabban/RosettaAbinitio) script to predict the fold of the desgined structure.
6. NOTE: Some backbones are very difficult to simulate their folds (*Abinitio*) before/after design becuase of a small energy gap between its fold and alternative folds using the same sequence, it is thus advised that you choose many structures and simulate them until you find a hit.
7. Generate fragments from the [Robetta server](http://robetta.org/) using this command:

`python3 Fragments.py FILENAME.pdb USERNAME`

This will generate and download the nessesary files to run an *Abinitio* simulation: FILENAME.pdb, structure.fasta, frags.200.3mers, frags.200.9mers, pre.psipred.ss2, it will also calculate the quality of the fragments and plot them in the file plot_frag.pdf.
