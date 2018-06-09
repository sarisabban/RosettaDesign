#!/usr/bin/python3
import os
import sys
import Bio.PDB
from Bio import pairwise2
from pyrosetta import *
from pyrosetta.toolbox import *
init()

class RosettaDesign():
	'''
	This class preforms RosettaDesign either fixed backbone 
	design (fixbb) or flexible backbone design (flxbb).
	It is preferred to perform the design many times and 
	select the best (lowest) scoring structure.
	'''
	def __init__(self):
		pass

	def BLAST(self , filename1 , filename2):
		'''
		Performs a BLAST alignment between two sequences and prints
		the sequences as well as the percentage of sequence
		similarity
		'''
		seq1 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('filename1' , filename1) , aa_only = True)[0].get_sequence()
		seq2 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('filename2' , filename2) , aa_only = True)[0].get_sequence()
		alignment = pairwise2.align.globalxx(seq1 , seq2)
		total = alignment[0][4]
		similarity = alignment[0][2]
		percentage = (similarity * 100) / total
		print(seq1)
		print(seq2)
		print('Sequence Similarity: {}%'.format(percentage))

	def fixbb(self , filename , relax_iters , design_iters):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while maintaining a fixed backbone.
		Generates the structure.pdb file
		'''
		#A - Relax original structure
		pose = pose_from_pdb(filename)
		chain = pose.pdb_info().chain(1)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		Rscore_before = scorefxn(pose)
		Rpose_work = Pose()
		Rpose_lowest = Pose()
		Rscores = []
		Rscores.append(Rscore_before)
		for nstruct in range(relax_iters):
			Rpose_work.assign(pose)
			relax.apply(Rpose_work)
			Rscore_after = scorefxn(Rpose_work)
			Rscores.append(Rscore_after)
			if Rscore_after < Rscore_before:
				Rscore_before = Rscore_after
				Rpose_lowest.assign(Rpose_work)
			else:
				continue
		pose.assign(Rpose_lowest)
		RFinalScore = scorefxn(pose)
		#B - Perform fixbb RosettaDesign
		packtask = standard_packer_task(pose)
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , packtask)
		backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
		backrub.pivot_residues(pose)
		GMC = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover()
		GMC.set_mover(backrub)
		GMC.set_scorefxn(scorefxn)
		GMC.set_maxtrials(500)
		GMC.set_temperature(1.0)
		GMC.set_preapply(False)
		GMC.set_recover_low(True)
		mover = pyrosetta.rosetta.protocols.moves.SequenceMover()
		mover.add_mover(pack)
		mover.add_mover(GMC)#####<--- problem here not accepting not rejecting moves
		Dscore_before = 0
		Dpose_work = Pose()
		Dpose_lowest = Pose()
		Dscores = []
		Dscores.append(Dscore_before)
		for nstruct in range(design_iters):
			Dpose_work.assign(pose)
			mover.apply(Dpose_work)
			Dscore_after = scorefxn(Dpose_work)
			Dscores.append(Dscore_after)
			if Dscore_after < Dscore_before:
				Dscore_before = Dscore_after
				Dpose_lowest.assign(Dpose_work)
			else:
				continue
		pose.assign(Dpose_lowest)
		DFinalScore = scorefxn(pose)
		#C - Output Result
		pose.dump_pdb('structure.pdb')
		#D - Print report
		print('==================== Result Report ====================')
		print('Relax Scores:\n' , Rscores)
		print('Chosen Lowest Score:' , RFinalScore , '\n')
		print('Design Scores:\n' , Dscores)
		print('Chosen Lowest Score:' , DFinalScore , '\n')
		print('BLAST result, compairing the original structure to the designed structure:')
		RosettaDesign.BLAST(self , filename , 'structure.pdb')

	def flxbb(self , filename , relax_iters , design_iters):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while allowing for a flexible backbone.
		Generates the structure.pdb file
		'''

def main(protocol , filename):
	RD = RosettaDesign()
	if protocol == 'fixbb':
		RD.fixbb(filename , 50 , 100)	##### <-------- Requires Work
	elif protocol == 'flxbb':
		RD.flxbb(filename , 50 , 100)	##### <-------- Requires Work

if __name__ == '__main__':
	main(sys.argv[1] , sys.argv[2])