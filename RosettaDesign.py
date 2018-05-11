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

	def SASA(self , filename):
		'''
		Calculates the different layers 
		(Surface, Boundary, Core) of a structure according
		its SASA (solvent-accessible surface area)
		Returns three lists
		Surface amino acids = [0]
		Boundary amino acids = [1]
		Core amino acids = [2]
		'''

		#Standard script to setup biopython's DSSP to calculate SASA using the Wilke constants
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('temp' , filename)
		dssp = Bio.PDB.DSSP(structure[0] , filename , acc_array = 'Wilke')
		#Loop to get SASA for each amino acid
		lis = list()
		count = 0
		for x in dssp:
			if x[1]   == 'A' : sasa = 129 * (x[3])
			elif x[1] == 'V' : sasa = 174 * (x[3])
			elif x[1] == 'I' : sasa = 197 * (x[3])
			elif x[1] == 'L' : sasa = 201 * (x[3])
			elif x[1] == 'M' : sasa = 224 * (x[3])
			elif x[1] == 'P' : sasa = 159 * (x[3])
			elif x[1] == 'Y' : sasa = 263 * (x[3])
			elif x[1] == 'F' : sasa = 240 * (x[3])
			elif x[1] == 'W' : sasa = 285 * (x[3])
			elif x[1] == 'R' : sasa = 274 * (x[3])
			elif x[1] == 'C' : sasa = 167 * (x[3])
			elif x[1] == 'N' : sasa = 195 * (x[3])
			elif x[1] == 'Q' : sasa = 225 * (x[3])
			elif x[1] == 'E' : sasa = 223 * (x[3])
			elif x[1] == 'G' : sasa = 104 * (x[3])
			elif x[1] == 'H' : sasa = 224 * (x[3])
			elif x[1] == 'K' : sasa = 236 * (x[3])
			elif x[1] == 'S' : sasa = 155 * (x[3])
			elif x[1] == 'T' : sasa = 172 * (x[3])
			elif x[1] == 'D' : sasa = 193 * (x[3])
			lis.append((x[2] , sasa))
		#Label each amino acid depending on its SASA position according to the parameters highlighted in the paper by (Koga et.al., 2012 - PMID: 23135467). The parameters are as follows:
		#Surface:	Helix or Sheet: SASA => 60		Loop: SASA => 40
		#Boundry:	Helix or Sheet: 15 < SASA < 60		Loop: 25 < SASA < 40
		#Core:		Helix or Sheet: SASA =< 15		Loop: SASA =< 25	
		surface = list()
		boundary = list()
		core = list()
		count = 0
		for x , y in lis:
			count = count + 1
			if y <= 25 and (x == '-' or x == 'T' or x == 'S'):					#Loop (DSSP code is - or T or S)
				core.append(count)
			elif 25 < y < 40 and (x == '-' or x == 'T' or x == 'S'):				#Loop (DSSP code is - or T or S)
				boundary.append(count)
			elif y >= 40 and (x == '-' or x == 'T' or x == 'S'):					#Loop (DSSP code is - or T or S)
				surface.append(count)
			elif y <= 15 and (x == 'G' or x == 'H' or x == 'I'):					#Helix (DSSP code is G or H or I)
				core.append(count)
			elif 15 < y < 60 and (x == 'G' or x == 'H' or x == 'I'):				#Helix (DSSP code is G or H or I)
				boundary.append(count)
			elif y >= 60 and (x == 'G' or x == 'H' or x == 'I'):					#Helix (DSSP code is G or H or I)
				surface.append(count)
			elif y <= 15 and (x == 'B' or x == 'E'):						#Sheet (DSSP code is B or E)
				core.append(count)
			elif 15 < y < 60 and (x == 'B' or x == 'E'):						#Sheet (DSSP code is B or E)
				boundary.append(count)
			elif y >= 60 and (x == 'B' or x == 'E'):						#Sheet (DSSP code is B or E)
				surface.append(count)	
		return(surface , boundary , core)

	def BLAST(self , filename1 , filename2):
		seq1 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('filename1' , filename1) , aa_only = True)[0].get_sequence()
		seq2 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('filename2' , filename2) , aa_only = True)[0].get_sequence()
		alignment = pairwise2.align.globalxx(seq1 , seq2)
		total = alignment[0][4]
		similarity = alignment[0][2]
		percentage = (similarity * 100) / total
		print(seq1)
		print(seq2)
		print('Sequence Similarity: {}%'.format(percentage))

	def whole_fixbb(self , filename):
		'''
		Applies RosettaDesign to change the whole
		structure's amino acids (the whole structure all
		at once) while maintaining the same backbone
		Generates the structure.pdb file
		'''
		#A - Relax original structure
		pose = pose_from_pdb(filename)
		scorefxn = get_fa_scorefxn()
		score1_original_before_relax = scorefxn(pose)							#Measure score before relaxing
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		relax.apply(pose)
		score2_original_after_relax = scorefxn(pose)							#Measure score after relaxing
		#B - Preform RosettaDesign for whole structure
		for inter in range(3):
			task_pack = standard_packer_task(pose)
			pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task_pack)
			pack_mover.apply(pose)
			#C - Relax Pose
			relax.apply(pose)
		#D - Output Result
		score3_of_design_after_relax = scorefxn(pose)							#Measure score of designed pose
		pose.dump_pdb('structure.pdb')									#Export final pose into a .pdb structure file
		print('---------------------------------------------------------')
		print('Original Structure Score:\t{}'.format(score1_original_before_relax))
		print('Relaxed Original Score:\t{}'.format(score2_original_after_relax))
		print('Relaxed Design Score:\t{}'.format(score3_of_design_after_relax))
		RosettaDesign.BLAST(self , filename , 'structure.pdb')

	def layer_fixbb(self , filename):
		'''
		Applies RosettaDesign to change the whole
		structure's amino acids (one layer at a time) while
		maintaining the same backbone. Should be more
		efficient and faster than the previous Whole method
		Generates the structure.pdb file
		'''
		#A - Relax original structure
		pose = pose_from_pdb(filename)
		scorefxn = get_fa_scorefxn()
		score1_original_before_relax = scorefxn(pose)							#Measure score before relaxing
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		relax.apply(pose)
		score2_original_after_relax = scorefxn(pose)							#Measure score after relaxing
		#B - Preform RosettaDesign one layer at a time
		for inter in range(3):
			#1 - Get SASA Layers
			sasa = RosettaDesign.SASA(self , filename)
			surface = sasa[0]
			boundary = sasa[1]
			core = sasa[2]
			#2 - Perform RosettaDesign on each layer
			#Design core
			task_pack = standard_packer_task(pose)
			pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , task_pack)
			task_pack.temporarily_fix_everything()							#To prevent all amino acids from being designed
			for AA in core:
				coreAA = pose.residue(AA).name()
				if coreAA == 'CYS:disulfide':
					continue
				else:
					task_pack.temporarily_set_pack_residue(AA , True)			#To move only spesific amino acids
			pack_mover.apply(pose)
			#Design boundery
			task_pack = standard_packer_task(pose)
			pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , task_pack)
			task_pack.temporarily_fix_everything()							#To prevent all amino acids from being designed
			for AA in boundary:
				boundAA = pose.residue(AA).name()
				if boundAA == 'CYS:disulfide':
					continue
				else:
					task_pack.temporarily_set_pack_residue(AA , True)			#To move only spesific amino acids
			pack_mover.apply(pose)
			#Design surface
			task_pack = standard_packer_task(pose)
			pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn , task_pack)
			task_pack.temporarily_fix_everything()							#To prevent all amino acids from being designed
			for AA in surface:
				surfAA = pose.residue(AA).name()
				if surfAA == 'CYS:disulfide':
					continue
				else:
					task_pack.temporarily_set_pack_residue(AA , True)			#To move only spesific amino acids
			pack_mover.apply(pose)
			#3 - Relax Pose
			relax.apply(pose)
		#C - Output Result
		score3_of_design_after_relax = scorefxn(pose)							#Measure score of designed pose
		pose.dump_pdb('structure.pdb')									#Export final pose into a .pdb structure file
		print('---------------------------------------------------------')
		print('Original Structure Score:\t{}'.format(score1_original_before_relax))
		print('Relaxed Original Score:\t{}'.format(score2_original_after_relax))
		print('Relaxed Design Score:\t{}'.format(score3_of_design_after_relax))
		RosettaDesign.BLAST(self , filename , 'structure.pdb')

	def pack_fixbb(self , filename):
		'''
		Applies RosettaDesign to change the whole
		structure's amino acids (one layer at a time as
		well as designing towards an optimally packed core)
		while maintaining the same backbone
		Generates the structure.pdb file
		'''
		#A - Relax original structure
		pose = pose_from_pdb(filename)
		scorefxn = get_fa_scorefxn()
		score1_original_before_relax = scorefxn(pose)							#Measure score before relaxing
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		relax.apply(pose)
		score2_original_after_relax = scorefxn(pose)							#Measure score after relaxing
		#B - FastRelax protocol										#Uses Generic Monte Carlo with PackStat as a filter to direct FastRelax towards an optimally packed structure core
		chain = pose.pdb_info().chain(1)								#Identify chain
		layers = [2 , 1 , 0]										#Layer Identity from SASA Surface = [0] , Boundary = [1] , Core = [2]
		for identity in layers:										#Loop through each layer
			#1 - Setup the PackStat filter
			filters = rosetta.protocols.simple_filters.PackStatFilter()
			#2 - Identify The Layers
			sasa = RosettaDesign.SASA(self , filename)						#Re-calculate SASA every time because amino acid position can change from one layer to another during the design phase, therefore make sure to design the layer not the amino acid
			layer = sasa[identity]									#Changes every iteration to start with Core (sasa[2]) then Boundary (sasa[1]) then Surface (sasa[0])
			#3 - Generate the resfile								#Will generate a new Resfile for each layer (which is why it is deleted at the end of the loop)
			Resfile = open('Resfile.resfile' , 'w')
			Resfile.write('NATAA\n')
			Resfile.write('start\n')
			for line in layer:
				Resfile.write(str(line) + ' ' + chain + ' ALLAA\n')
			Resfile.close()
			#4 - Setup the FastRelax mover
			read = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('Resfile.resfile')	#Call the generated Resfile
			task = pyrosetta.rosetta.core.pack.task.TaskFactory()					#Setup the TaskFactory
			task.push_back(read)									#Add the Resfile to the TaskFactory
			movemap = MoveMap()									#Setup the MoveMap
			movemap.set_bb(False)									#Do not change the phi and psi BackBone angles
			movemap.set_chi(True)									#Change the chi Side Chain angle
			mover = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()			#Call the FastDesign mover
			mover.set_task_factory(task)								#Add the TaskFactory to it
			mover.set_movemap(movemap)								#Add the MoveMap to it
			mover.set_scorefxn(scorefxn)								#Add the Score Function to it
			mover.apply(pose)
			os.remove('Resfile.resfile')								#To keep working directory clean, and to make sure each Resfile has the info for each layer only and they do not get mixed and appended together in one Resfile
			#5 - Relax pose
			relax.apply(pose)
		#C - Output result
		score3_of_design_after_relax = scorefxn(pose)							#Measure score of designed pose
		pose.dump_pdb('structure.pdb')									#Export final pose into a .pdb structure file
		print('---------------------------------------------------------')
		print('Original Structure Score:\t{}'.format(score1_original_before_relax))
		print('Relaxed Original Score:\t{}'.format(score2_original_after_relax))
		print('Relaxed Design Score:\t{}'.format(score3_of_design_after_relax))
		RosettaDesign.BLAST(self , filename , 'structure.pdb')

	def pack_flxbb(self , filename):
		'''
		Applies RosettaDesign to change the whole
		structure's amino acids (one layer at a time as
		well as designing towards an optimally packed core)
		while maintaining the same backbone
		Generates the structure.pdb file
		'''
		#A - Relax original structure
		pose = pose_from_pdb(filename)
		scorefxn = get_fa_scorefxn()
		score1_original_before_relax = scorefxn(pose)							#Measure score before relaxing
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		relax.apply(pose)
		score2_original_after_relax = scorefxn(pose)							#Measure score after relaxing
		#B - FastRelax protocol										#Uses Generic Monte Carlo with PackStat as a filter to direct FastRelax towards an optimally packed structure core
		chain = pose.pdb_info().chain(1)								#Identify chain
		layers = [2 , 1 , 0]										#Layer Identity from SASA Surface = [0] , Boundary = [1] , Core = [2]
		for identity in layers:										#Loop through each layer
			#1 - Setup the PackStat filter
			filters = rosetta.protocols.simple_filters.PackStatFilter()
			#2 - Identify The Layers
			sasa = RosettaDesign.SASA(self , filename)						#Re-calculate SASA every time because amino acid position can change from one layer to another during the design phase, therefore make sure to design the layer not the amino acid
			layer = sasa[identity]									#Changes every iteration to start with Core (sasa[2]) then Boundary (sasa[1]) then Surface (sasa[0])
			#3 - Generate the resfile								#Will generate a new Resfile for each layer (which is why it is deleted at the end of the loop)
			Resfile = open('Resfile.resfile' , 'w')
			Resfile.write('NATAA\n')
			Resfile.write('start\n')
			for line in layer:
				Resfile.write(str(line) + ' ' + chain + ' ALLAA\n')
			Resfile.close()
			#4 - Setup the FastRelax mover
			read = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('Resfile.resfile')	#Call the generated Resfile
			task = pyrosetta.rosetta.core.pack.task.TaskFactory()					#Setup the TaskFactory
			task.push_back(read)									#Add the Resfile to the TaskFactory
			movemap = MoveMap()									#Setup the MoveMap
			movemap.set_bb(True)									#Do not change the phi and psi BackBone angles
			movemap.set_chi(True)									#Change the chi Side Chain angle
			mover = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()			#Call the FastDesign mover
			mover.set_task_factory(task)								#Add the TaskFactory to it
			mover.set_movemap(movemap)								#Add the MoveMap to it
			mover.set_scorefxn(scorefxn)								#Add the Score Function to it
			mover.apply(pose)
			os.remove('Resfile.resfile')								#To keep working directory clean, and to make sure each Resfile has the info for each layer only and they do not get mixed and appended together in one Resfile
			#5 - Relax pose
			relax.apply(pose)
		#C - Output result
		score3_of_design_after_relax = scorefxn(pose)							#Measure score of designed pose
		pose.dump_pdb('structure.pdb')									#Export final pose into a .pdb structure file
		print('---------------------------------------------------------')
		print('Original Structure Score:\t{}'.format(score1_original_before_relax))
		print('Relaxed Original Score:\t{}'.format(score2_original_after_relax))
		print('Relaxed Design Score:\t{}'.format(score3_of_design_after_relax))
		RosettaDesign.BLAST(self , filename , 'structure.pdb')

def main():
	protocol = sys.argv[1]
	filename = sys.argv[2]

	RD = RosettaDesign()

	if protocol == 'fixbb':
		#RD.whole_fixbb(filename)
		#RD.layer_fixbb(filename)
		RD.pack_fixbb(filename)
		#RD.pack_flxbb(filename)	

	elif protocol == 'flxbb':
		#RD.whole_fixbb(filename)
		#RD.layer_fixbb(filename)
		#RD.pack_fixbb(filename)
		RD.pack_flxbb(filename)

if __name__ == '__main__':
	main()
