#!/usr/bin/python3

import numpy , Bio.PDB , os , sys

from pyrosetta import *
from pyrosetta.toolbox import *
from rosetta.protocols.relax import *
from rosetta.protocols.moves import *
from rosetta.protocols.simple_moves import *
from rosetta.core.pack.task import TaskFactory
init()

def Relax(pose):
	''' Relaxes a structure '''
	''' Updates the original pose with the relaxed pose '''
	scorefxn = get_fa_scorefxn()
	relax = ClassicRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(pose)

def SASA(pose):
	''' Calculates the different layers (Surface, Boundery, Core) of a structure according its SASA (solvent-accessible surface area) '''
	''' Returns three lists Surface amino acids = [0] , Boundery amino acids = [1] , Core amino acids = [2] '''
	#Temporary generate a .pdb file of the pose to isolate the layers since it is not yet possible to do that using a Rosetta pose, this temporary .pdb file will be deleted after the layers are found
	pose.dump_pdb('ToDesign.pdb')
	#Standard script to setup biopython's DSSP to calculate SASA using Wilke constants
	p = Bio.PDB.PDBParser()
	structure = p.get_structure('X' , 'ToDesign.pdb')
	model = structure[0]
	dssp = Bio.PDB.DSSP(model , 'ToDesign.pdb' , acc_array='Wilke')
	#Loop to get SASA for each amino acid
	lis = list()
	count = 0
	for x in dssp:
		if x[1]=='A' : sasa=129*(x[3])
		elif x[1]=='V' : sasa=174*(x[3])
		elif x[1]=='I' : sasa=197*(x[3])
		elif x[1]=='L' : sasa=201*(x[3])
		elif x[1]=='M' : sasa=224*(x[3])
		elif x[1]=='P' : sasa=159*(x[3])
		elif x[1]=='Y' : sasa=263*(x[3])
		elif x[1]=='F' : sasa=240*(x[3])
		elif x[1]=='W' : sasa=285*(x[3])
		elif x[1]=='R' : sasa=274*(x[3])
		elif x[1]=='C' : sasa=167*(x[3])
		elif x[1]=='N' : sasa=195*(x[3])
		elif x[1]=='Q' : sasa=225*(x[3])
		elif x[1]=='E' : sasa=223*(x[3])
		elif x[1]=='G' : sasa=104*(x[3])
		elif x[1]=='H' : sasa=224*(x[3])
		elif x[1]=='K' : sasa=236*(x[3])
		elif x[1]=='S' : sasa=155*(x[3])
		elif x[1]=='T' : sasa=172*(x[3])
		elif x[1]=='D' : sasa=193*(x[3])
		lis.append((x[2] , sasa))
	#Label each amino acid depending on its SASA position according to the parameters highlighted in the paper by (Koga et.al., 2012 - PMID: 23135467). The parameters are as follows:
	#Surface:	Helix or Sheet: SASA => 60		Loop: SASA => 40
	#Boundery:	Helix or Sheet: 15 < SASA < 60		Loop: 25 < SASA < 40
	#Core:		Helix or Sheet: SASA =< 15		Loop: SASA =< 25	
	surface = list()
	boundery = list()
	core = list()
	count = 0
	for x , y in lis:
		count = count + 1
		if y <= 25 and (x == '-' or x == 'T' or x == 'S'):			#Loop (DSSP code is - or T or S)
			core.append(count)
		elif 25 < y < 40 and (x == '-' or x == 'T' or x == 'S'):		#Loop (DSSP code is - or T or S)
			boundery.append(count)
		elif y >= 40 and (x == '-' or x == 'T' or x == 'S'):			#Loop (DSSP code is - or T or S)
			surface.append(count)
		elif y <= 15 and (x == 'G' or x == 'H' or x == 'I'):			#Helix (DSSP code is G or H or I)
			core.append(count)
		elif 15 < y < 60 and (x == 'G' or x == 'H' or x == 'I'):		#Helix (DSSP code is G or H or I)
			boundery.append(count)
		elif y >= 60 and (x == 'G' or x == 'H' or x == 'I'):			#Helix (DSSP code is G or H or I)
			surface.append(count)
		elif y <= 15 and (x == 'B' or x == 'E'):				#Sheet (DSSP code is B or E)
			core.append(count)
		elif 15 < y < 60 and (x == 'B' or x == 'E'):				#Sheet (DSSP code is B or E)
			boundery.append(count)
		elif y >= 60 and (x == 'B' or x == 'E'):				#Sheet (DSSP code is B or E)
			surface.append(count)	
	os.remove('ToDesign.pdb')							#Keep working directory clean
	return(surface , boundery , core)

def Design_Whole(pose):
	''' Applies RosettaDesign to change the whole structure's amino acids (the whole structure all at once) while maintaining the same backbone '''
	''' Generates the Designed.pdb file '''
	#1 - Relax original structure
	scorefxn = get_fa_scorefxn()							#Call the score function
	score1_original_before_relax = scorefxn(pose)					#Measure score before relaxing
	Relax(pose)									#Relax structure
	score2_original_after_relax = scorefxn(pose)					#Measure score after relaxing
	#2 - Preform RosettaDesign for whole structure
	for inter in range(3):
		task_pack = standard_packer_task(pose)
		pack_mover = PackRotamersMover(scorefxn, task_pack)
		pack_mover.apply(pose)
		#3 - Relax pose
		Relax(pose)
	score3_of_design_after_relax = scorefxn(pose)					#Measure score of designed pose
	pose.dump_pdb('Designed.pdb')							#Export final pose into a .pdb structure file
	print(score1_original_before_relax)
	print(score2_original_after_relax)
	print(score3_of_design_after_relax)

def Design_Layer(pose):
	''' Applies RosettaDesign to change the whole structure's amino acids (one layer at a time) while maintaining the same backbone. Should be more efficient and faster than the previous method (Design_Full) '''
	''' Generates the Designed.pdb file '''
	#Relax original structure
	scorefxn = get_fa_scorefxn()							#Call the score function
	score1_original_before_relax = scorefxn(pose)					#Measure score before relaxing
	Relax(pose)									#Relax structure
	score2_original_after_relax = scorefxn(pose)					#Measure score after relaxing
	#Preform RosettaDesign one layer at a time
	for inter in range(3):
		#1 - Get SASA Layers
		sasa = SASA(pose)
		surface = sasa[0]
		boundery = sasa[1]
		core = sasa[2]
		#2 - Preform RosettaDesign on each layer
		#Design core
		task_pack = standard_packer_task(pose)
		pack_mover = PackRotamersMover(scorefxn , task_pack)
		task_pack.temporarily_fix_everything()					#To prevent all amino acids from being designed
		for AA in core:
			coreAA = pose.residue(AA).name()
			if coreAA == 'CYS:disulfide':
				continue
			else:
				task_pack.temporarily_set_pack_residue(AA , True)	#To move only spesific amino acids
		pack_mover.apply(pose)
		#Design boundery
		task_pack = standard_packer_task(pose)
		pack_mover = PackRotamersMover(scorefxn , task_pack)
		task_pack.temporarily_fix_everything()					#To prevent all amino acids from being designed
		for AA in boundery:
			boundAA = pose.residue(AA).name()
			if boundAA == 'CYS:disulfide':
				continue
			else:
				task_pack.temporarily_set_pack_residue(AA , True)	#To move only spesific amino acids
		pack_mover.apply(pose)
		#Design surface
		task_pack = standard_packer_task(pose)
		pack_mover = PackRotamersMover(scorefxn , task_pack)
		task_pack.temporarily_fix_everything()					#To prevent all amino acids from being designed
		for AA in surface:
			surfAA = pose.residue(AA).name()
			if surfAA == 'CYS:disulfide':
				continue
			else:
				task_pack.temporarily_set_pack_residue(AA , True)	#To move only spesific amino acids
		pack_mover.apply(pose)
		#3 - Relax pose
		Relax(pose)
	score3_of_design_after_relax = scorefxn(pose)					#Measure score of designed pose
	pose.dump_pdb('Designed.pdb')							#Export final pose into a .pdb structure file
	print(score1_original_before_relax)
	print(score2_original_after_relax)
	print(score3_of_design_after_relax)
#------------------------------------------------------------------------------------------------------------------------------------
pose = pose_from_pdb(sys.argv[1])
#Design_Whole(pose)
Design_Layer(pose)
