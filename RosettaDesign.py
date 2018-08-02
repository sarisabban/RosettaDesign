#!/usr/bin/python3

import os
import sys
import Bio.PDB
from Bio import pairwise2
from pyrosetta import *
from pyrosetta.toolbox import *
init()

def BUH(filename):
	'''
	Prints out all the residues form a structure
	with unsatisfied hydrogen bonds.
	'''
	import json
	pose = pose_from_pdb(filename)
	scorefxn = get_fa_scorefxn()
	scorefxn(pose)
	bupc = pyrosetta.rosetta.protocols.simple_pose_metric_calculators.BuriedUnsatisfiedPolarsCalculator('default', 'default')
	residues = bupc.get('residue_bur_unsat_polars', pose)
	buhs_for_each_res = json.loads(bupc.get('residue_bur_unsat_polars', pose))
	for i in range(pose.size()):
		if buhs_for_each_res[i] > 0:
			print('{0}{1}: {2}'.format(pose.pdb_info().chain(i + 1), pose.pdb_info().number(i + 1), buhs_for_each_res[i]))

class RosettaDesign():
	'''
	This class preforms RosettaDesign either fixed backbone
	design (fixbb) or flexible backbone design (flxbb).
	It is preferred to perform the design many times and 
	select the best (lowest) scoring structure.
	'''
	def __init__(self):
		pass

	def BLAST(self, filename1, filename2):
		'''
		Performs a BLAST alignment between two sequences and prints
		the sequences as well as the percentage of sequence
		similarity
		'''
		seq1 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET=True).get_structure('filename1', filename1), aa_only=True)[0].get_sequence()
		seq2 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET=True).get_structure('filename2', filename2), aa_only=True)[0].get_sequence()
		alignment = pairwise2.align.globalxx(seq1, seq2)
		total = alignment[0][4]
		similarity = alignment[0][2]
		percentage = (similarity*100)/total
		print(seq1)
		print(seq2)
		print('Sequence Similarity: {}%'.format(round(percentage, 3)))

	def fixbb(self, filename, relax_iters, design_iters):
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
		pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask)
		Dscore_before = 0
		Dpose_work = Pose()
		Dpose_lowest = Pose()
		Dscores = []
		Dscores.append(Dscore_before)
		for nstruct in range(design_iters):
			Dpose_work.assign(pose)
			pack.apply(Dpose_work)
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
		pose.dump_pdb('fixbb.pdb')
		#D - Print report
		print('==================== Result Report ====================')
		print('Relax Scores:\n', Rscores)
		print('Chosen Lowest Score:', RFinalScore, '\n')
		print('Design Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')
		print('BLAST result, comparing the original structure to the designed structure:')
		RosettaDesign.BLAST(self, filename, 'fixbb.pdb')

	def flxbb(self, filename, relax_iters, design_iters):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while allowing for a flexible backbone.
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
		#B - Perform flxbb RosettaDesign
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		mover = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		mover.set_task_factory(task)
		mover.set_movemap(movemap)
		mover.set_scorefxn(scorefxn)
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
		pose.dump_pdb('flxbb.pdb')
		#D - Print report
		print('==================== Result Report ====================')
		print('Relax Scores:\n', Rscores)
		print('Chosen Lowest Score:', RFinalScore, '\n')
		print('Design Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')
		print('BLAST result, comparing the original structure to the designed structure:')
		RosettaDesign.BLAST(self, filename, 'flxbb.pdb')

	def BDR(self, filename, refine_iters):
		#A - Generate constraints file
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure('{}'.format(filename), filename)
		length = len(structure[0]['A'])
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure, aa_only=False)
		model = Type
		chain = model[0]
		CST = []
		CST.append(0.0)
		for aa in range(1, length+1):
			try:
				residue1 = chain[0]
				residue2 = chain[aa]
				atom1 = residue1['CA']
				atom2 = residue2['CA']
				CST.append(atom1-atom2)
			except:
				pass
		atom = 1
		for cst in CST:
			line = 'AtomPair CA 1 CA '+str(atom)+' GAUSSIANFUNC '+str(cst)+' 1.0\n'
			thefile = open('structure.constraints', 'a')
			thefile.write(line)
			thefile.close()
			atom += 1
		#B - Generate blueprint file (remodeling only large loops)
		dssp = Bio.PDB.DSSP(structure[0], filename)
		SS = []
		SEQ = []
		for ss in dssp:
			if ss[2] == 'G' or ss[2] == 'H' or ss[2] == 'I':
				rename = 'HX'
			elif ss[2] == 'B' or ss[2] == 'E':
				rename = 'EX'
			else:
				rename = 'LX'
			SS.append(rename)
			SEQ.append(ss[1])
		buf = []
		items = []
		l_seen = 0
		for count, (ss, aa) in enumerate(zip(SS, SEQ), 1):
			buf.append((count, aa, ss))
			if 'LX' in {ss, aa}:
				l_seen += 1
				if l_seen >= 3:
					for count, aa, ss in buf:
						line = [str(count), aa, ss, '.' if ss in {'HX', 'EX'} else 'R']
						line = ' '.join(line)
						items.append(line)
					buf.clear()
			else:
				l_seen = 0
				for count, aa, ss in buf:
					line = [str(count), aa, ss, '.']
					line = ' '.join(line)
					items.append(line)
				buf.clear()
		if int(items[-1].split()[0]) != count:
			line = [str(count), aa, ss, '.']
			line = ' '.join(line)
			items.append(line)
		blueprint = open('structure.blueprint', 'a')
		for line in items:
			blueprint.write(line+'\n')
		blueprint.close()
		#C - Run BluePrint mover
		pose = pose_from_pdb(filename)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		secstr = pyrosetta.rosetta.protocols.fldsgn.potentials.SetSecStructEnergies(scorefxn,'structure.blueprint', True)
		secstr.apply(pose)
		BDR = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
		BDR.num_fragpick(200)
		BDR.use_fullmer(True)
		BDR.use_sequence_bias(False)
		BDR.max_linear_chainbreak(0.07)
		BDR.ss_from_blueprint(True)
		BDR.dump_pdb_when_fail('')
		BDR.set_constraints_NtoC(-1.0)
		BDR.use_abego_bias(True)
		#BDR.set_constraint_file('structure.constraints')
		BDR.set_blueprint('structure.blueprint')
		Dscore_before = 0
		Dpose_work = Pose()
		Dpose_lowest = Pose()
		Dscores = []
		Dscores.append(Dscore_before)
		for nstruct in range(refine_iters):
			Dpose_work.assign(pose)
			BDR.apply(Dpose_work)
			relax.apply(Dpose_work)
			Dscore_after = scorefxn(Dpose_work)
			Dscores.append(Dscore_after)
			if Dscore_after < Dscore_before:
				Dscore_before = Dscore_after
				Dpose_lowest.assign(Dpose_work)
			else:
				continue
		pose.assign(Dpose_lowest)
		DFinalScore = scorefxn(pose)
		#D - Output Result
		pose.dump_pdb('remodel.pdb')
		os.remove('structure.constraints')
		os.remove('structure.blueprint')
		#E - Print report
		print('==================== Result Report ====================')
		print('Design Scores:\n', Dscores)
		print('Chosen Lowest Score:', DFinalScore, '\n')

	def Refine(self, filename, refine_iters):
		os.system('cp {} temp.pdb'.format(filename))
		Mutate = [1]
		while Mutate != []:
			inputfile = 'temp.pdb'
			parser = Bio.PDB.PDBParser()
			structure = parser.get_structure('{}'.format(inputfile), inputfile)
			dssp = Bio.PDB.DSSP(structure[0], inputfile, acc_array='Wilke')
			sasalist = []
			for x in dssp:
				if x[1] == 'A':
					sasa = 129*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'V':
					sasa = 174*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'I':
					sasa = 197*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'L':
					sasa = 201*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'M':
					sasa = 224*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'P':
					sasa = 159*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'Y':
					sasa = 263*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'F':
					sasa = 240*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'W':
					sasa = 285*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'R':
					sasa = 274*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'N':
					sasa = 195*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'C':
					sasa = 167*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'Q':
					sasa = 225*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'E':
					sasa = 223*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'G':
					sasa = 104*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'H':
					sasa = 224*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'K':
					sasa = 236*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'S':
					sasa = 155*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'T':
					sasa = 172*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				elif x[1] == 'D':
					sasa = 193*(x[3])
					if sasa <= 25:
						sasa = 'C'
					elif 25 < sasa < 40:
						sasa = 'B'
					elif sasa >= 40:
						sasa = 'S'
				if x[2] == 'G' or x[2] == 'H' or x[2] == 'I':
					ss = 'H'
				elif x[2] == 'B' or x[2] == 'E':
					ss = 'S'
				elif x[2] == 'S' or x[2] == 'T' or x[2] == '-':
					ss = 'L'
				sasalist.append((x[0], x[1], ss, sasa))
			Resids = []
			SecStr = []
			SASAps = []
			MutPos = []
			Mutate = []
			for n, r, s, a in sasalist:
				if a == 'S' and s == 'L' and (	   r == 'P' or r == 'G' 
								or r == 'N' or r == 'Q'
								or r == 'S' or r == 'T'
								or r == 'D' or r == 'E'
								or r == 'R' or r == 'K'
								or r == 'H'):
					MutPos.append(' ')
				elif a=='B' and s=='L' and (	   r == 'A' or r == 'V'
								or r == 'I' or r == 'L'
								or r == 'F' or r == 'Y'
								or r == 'W' or r == 'G'
								or r == 'N' or r == 'Q'
								or r == 'S' or r == 'T'
								or r == 'P' or r == 'D'
								or r == 'E' or r == 'K'
								or r == 'R'):
					MutPos.append(' ')
				elif a=='C' and s=='L' and (	   r == 'A' or r == 'V'
								or r == 'I' or r == 'L'
								or r == 'P' or r == 'F'
								or r == 'W' or r == 'M'):
					MutPos.append(' ')
				elif a=='S' and s=='H' and (	   r == 'Q' or r == 'E'
								or r == 'K' or r == 'H'):
					MutPos.append(' ')
				elif a=='B' and s=='H' and (	   r == 'A' or r == 'V'
								or r == 'I' or r == 'L'
								or r == 'W' or r == 'Q'
								or r == 'E' or r == 'K'
								or r == 'F' or r == 'M'):
					MutPos.append(' ')
				elif a=='C' and s=='H' and (	   r == 'A' or r == 'V'
								or r == 'I' or r == 'L'
								or r == 'F' or r == 'W'
								or r == 'M'):
					MutPos.append(' ')
				elif a=='S' and s=='S' and (	   r == 'Q' or r == 'T'
								or r == 'Y'):
					MutPos.append(' ')
				elif a=='B' and s=='S' and (	   r == 'A' or r == 'V'
								or r == 'I' or r == 'L'
								or r == 'F' or r == 'Y'
								or r == 'W' or r == 'Q'
								or r == 'T' or r == 'M'):
					MutPos.append(' ')
				elif a=='C' and s=='S' and (	   r == 'A' or r == 'V'
								or r == 'I' or r == 'L'
								or r == 'F' or r == 'W'
								or r == 'M'):
					MutPos.append(' ')
				else:
					MutPos.append('*')
					Mutate.append((n, r, s, a))		
				Resids.append(r)
				SASAps.append(a)
				SecStr.append(s)
			Resids=''.join(Resids)
			SASAps=''.join(SASAps)
			MutPos=''.join(MutPos)
			SecStr=''.join(SecStr)
			print('{}\n{}\n{}\n{}'.format(Resids, SecStr, SASAps, MutPos))
			pose = pose_from_pdb(inputfile)
			scorefxn = get_fa_scorefxn()
			relax = pyrosetta.rosetta.protocols.relax.FastRelax()
			relax.set_scorefxn(scorefxn)
			ideal = pyrosetta.rosetta.protocols.idealize.IdealizeMover()
			resfile = open('structure.res', 'a')
			resfile.write('NATRO\nSTART\n')
			for n, r, a, s in Mutate:
				if s == 'S' and a == 'L':
					line = '{} A PIKAA PGNQSTDERKH\n'.format(n)
					resfile.write(line)
				elif s == 'S' and a == 'H':
					line = '{} A PIKAA QEKH\n'.format(n)
					resfile.write(line)
				elif s == 'S' and a == 'S':
					line = '{} A PIKAA QTY\n'.format(n)
					resfile.write(line)
				elif s == 'B' and a == 'L':
					line = '{} A PIKAA AVILFYWGNQSTPDEKR\n'.format(n)
					resfile.write(line)
				elif s == 'B' and a == 'H':
					line = '{} A PIKAA AVILWQEKFM\n'.format(n)
					resfile.write(line)
				elif s == 'B' and a == 'S':
					line = '{} A PIKAA AVILFYWQTM\n'.format(n)
					resfile.write(line)
				elif s == 'C' and a == 'L':
					line = '{} A PIKAA AVILPFWM\n'.format(n)
					resfile.write(line)
				elif s == 'C' and a == 'H':
					line = '{} A PIKAA AVILFWM\n'.format(n)
					resfile.write(line)
				elif s == 'C' and a == 'S':
					line = '{} A PIKAA AVILFWM\n'.format(n)
					resfile.write(line)
			resfile.close()
			pack = standard_packer_task(pose)
			pack.temporarily_fix_everything()
			pyrosetta.rosetta.core.pack.task.parse_resfile(pose, pack, 'structure.res')
			for n, r, s, a in Mutate:
				x = pose.residue(n).name()
				if x == 'CYS:disulphide':
					continue
				else:
					pack.temporarily_set_pack_residue(n, True) 
			print(pack)
			pack = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, pack)
			Dscore_before = 0
			Dpose_work = Pose()
			Dpose_lowest = Pose()
			Dscores = []
			Dscores.append(Dscore_before)
			for nstruct in range(refine_iters):
				Dpose_work.assign(pose)
				pack.apply(Dpose_work)
				ideal.apply(Dpose_work)
				relax.apply(Dpose_work)
				Dscore_after = scorefxn(Dpose_work)
				Dscores.append(Dscore_after)
				if Dscore_after < Dscore_before:
					Dscore_before = Dscore_after
					Dpose_lowest.assign(Dpose_work)
				else:
					continue
			pose.assign(Dpose_lowest)
			DFinalScore = scorefxn(pose)
			os.remove('structure.res')
			os.remove('temp.pdb')
			pose.dump_pdb('temp.pdb')
		pose.dump_pdb('structure.pdb')
		os.remove('temp.pdb')

class MCRosettaDesign():
	'''
	This class preforms RosettaDesign either fixed backbone 
	design (fixbb) or flexible backbone design (flxbb) using
	the Monte Carlo method, generating many designed structures
	thus it is best to select the lowest scoring structure.
	'''
	def __init__(self):
		pass

	def fixbb(self, filename, kT, cycles, jobs, job_output):
		'''
		Performs fixed backbone RosettaDesign using the
		Monte Carlo method using the following sequence:
		1. Relax
		2. Fixed backbone design (by SASA layers)
		'''
		# Generate resfile
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
		sasalist = []
		for x in dssp:
			if x[1] == 'A':
				sasa = 129*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'V':
				sasa = 174*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'I':
				sasa = 197*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'L':
				sasa = 201*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'M':
				sasa = 224*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'P':
				sasa = 159*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'Y':
				sasa = 263*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'F':
				sasa = 240*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'W':
				sasa = 285*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'R':
				sasa = 274*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'N':
				sasa = 195*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'C':
				sasa = 167*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'Q':
				sasa = 225*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'E':
				sasa = 223*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'G':
				sasa = 104*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'H':
				sasa = 224*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'K':
				sasa = 236*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'S':
				sasa = 155*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'T':
				sasa = 172*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'D':
				sasa = 193*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			if x[2] == 'G' or x[2] == 'H' or x[2] == 'I':
				ss = 'H'
			elif x[2] == 'B' or x[2] == 'E':
				ss = 'S'
			elif x[2] == 'S' or x[2] == 'T' or x[2] == '-':
				ss = 'L'
			sasalist.append((x[0], x[1], ss, sasa))
		resfile = open('resfile', 'a')
		resfile.write('NATRO\nSTART\n')
		for n, r, a, s in sasalist:
			if s == 'S' and a == 'L':
				line = '{} A PIKAA PGNQSTDERKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'H':
				line = '{} A PIKAA QEKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'S':
				line = '{} A PIKAA QTY\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'L':
				line = '{} A PIKAA AVILFYWGNQSTPDEKR\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'H':
				line = '{} A PIKAA AVILWQEKFM\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'S':
				line = '{} A PIKAA AVILFYWQTM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'L':
				line = '{} A PIKAA AVILPFWM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'H':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'S':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
		resfile.close()
		# RosettaDesign: Relax Fixbb, Relax
		pose = pose_from_pdb(filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxnBUH = get_fa_scorefxn()
		scorefxnBUH.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.buried_unsatisfied_penalty, 1.0)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		packtask = standard_packer_task(pose)
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose, packtask, 'resfile')
		fixbb = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxnBUH, packtask)
		sequence = SequenceMover()
		sequence.add_mover(relax)
		sequence.add_mover(fixbb)
		sequence.add_mover(relax)
		mc = MonteCarlo(pose, scorefxn, kT)
		trial = TrialMover(sequence, mc)
		RosettaDesign = RepeatMover(trial, cycles)
		job = PyJobDistributor(job_output, jobs, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			mc.reset(pose)
			RosettaDesign.apply(pose)
			mc.recover_low(pose)
			job.output_decoy(pose)
		os.remove('resfile')

	def flxbb(self, filename, kT, cycles, jobs, job_output):
		'''
		Performs flexible backbone RosettaDesign using the
		Monte Carlo method using the following sequence:
		1. Relax
		2. BluePrintBDR loop remodeling
		3. Flexible backbone design (by SASA layers)
		4. Idealise
		5. Relax
		'''
		# Generate blueprint file
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename)
		SS = []
		SEQ = []
		for ss in dssp:
			if ss[2] == 'G' or ss[2] == 'H' or ss[2] == 'I':
				rename = 'HX'
			elif ss[2] == 'B' or ss[2] == 'E':
				rename = 'EX'
			else:
				rename = 'LX'
			SS.append(rename)
			SEQ.append(ss[1])
		buf = []
		items = []
		l_seen = 0
		for count, (ss, aa) in enumerate(zip(SS, SEQ), 1):
			buf.append((count, aa, ss))
			if 'LX' in {ss, aa}:
				l_seen += 1
				if l_seen >= 3:
					for count, aa, ss in buf:
						line = [str(count), aa, ss, '.' if ss in {'HX', 'EX'} else 'R']
						line = ' '.join(line)
						items.append(line)
					buf.clear()
			else:
				l_seen = 0
				for count, aa, ss in buf:
					line = [str(count), aa, ss, '.']
					line = ' '.join(line)
					items.append(line)
				buf.clear()
		if int(items[-1].split()[0]) != count:
			line = [str(count), aa, ss, '.']
			line = ' '.join(line)
			items.append(line)
		blueprint = open('blueprint', 'a')
		for line in items:
			blueprint.write(line+'\n')
		blueprint.close()
		# Generate resfile
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
		sasalist = []
		for x in dssp:
			if x[1] == 'A':
				sasa = 129*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'V':
				sasa = 174*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'I':
				sasa = 197*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'L':
				sasa = 201*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'M':
				sasa = 224*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'P':
				sasa = 159*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'Y':
				sasa = 263*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'F':
				sasa = 240*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'W':
				sasa = 285*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'R':
				sasa = 274*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'N':
				sasa = 195*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'C':
				sasa = 167*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'Q':
				sasa = 225*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'E':
				sasa = 223*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'G':
				sasa = 104*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'H':
				sasa = 224*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'K':
				sasa = 236*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'S':
				sasa = 155*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'T':
				sasa = 172*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			elif x[1] == 'D':
				sasa = 193*(x[3])
				if sasa <= 25:
					sasa = 'C'
				elif 25 < sasa < 40:
					sasa = 'B'
				elif sasa >= 40:
					sasa = 'S'
			if x[2] == 'G' or x[2] == 'H' or x[2] == 'I':
				ss = 'H'
			elif x[2] == 'B' or x[2] == 'E':
				ss = 'S'
			elif x[2] == 'S' or x[2] == 'T' or x[2] == '-':
				ss = 'L'
			sasalist.append((x[0], x[1], ss, sasa))
		resfile = open('resfile', 'a')
		resfile.write('NATRO\nSTART\n')
		for n, r, a, s in sasalist:
			if s == 'S' and a == 'L':
				line = '{} A PIKAA PGNQSTDERKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'H':
				line = '{} A PIKAA QEKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'S':
				line = '{} A PIKAA QTY\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'L':
				line = '{} A PIKAA AVILFYWGNQSTPDEKR\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'H':
				line = '{} A PIKAA AVILWQEKFM\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'S':
				line = '{} A PIKAA AVILFYWQTM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'L':
				line = '{} A PIKAA AVILPFWM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'H':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'S':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
		resfile.close()
		# RosettaDesign: Relax, BluePrintBDR, Flxbb, Idealize, Relax
		pose = pose_from_pdb(filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxnBUH = get_fa_scorefxn()
		scorefxnBUH.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.buried_unsatisfied_penalty, 1.0)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		BDR = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
		BDR.num_fragpick(200)
		BDR.use_fullmer(True)
		BDR.use_sequence_bias(False)
		BDR.max_linear_chainbreak(0.07)
		BDR.ss_from_blueprint(True)
		BDR.dump_pdb_when_fail('')
		BDR.set_constraints_NtoC(-1.0)
		BDR.use_abego_bias(True)
		BDR.set_blueprint('blueprint')
		resfile = rosetta.core.pack.task.operation.ReadResfile('resfile')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		flxbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		flxbb.set_task_factory(task)
		flxbb.set_movemap(movemap)
		flxbb.set_scorefxn(scorefxnBUH)
		ideal = pyrosetta.rosetta.protocols.idealize.IdealizeMover()
		sequence = SequenceMover()
		sequence.add_mover(relax)
		sequence.add_mover(BDR)
		sequence.add_mover(flxbb)
		sequence.add_mover(ideal)
		sequence.add_mover(relax)
		mc = MonteCarlo(pose, scorefxn, kT)
		trial = TrialMover(sequence, mc)
		RosettaDesign = RepeatMover(trial, cycles)
		job = PyJobDistributor(job_output, jobs, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			mc.reset(pose)
			RosettaDesign.apply(pose)
			mc.recover_low(pose)
			job.output_decoy(pose)
		os.remove('blueprint')
		os.remove('resfile')

def Protocol(protocol, filename):
	RD = RosettaDesign()
	if protocol == 'fixbb':
		RD.BDR(filename, 200)
		RD.fixbb('remodel.pdb', 50, 100)
		RD.Refine('fixbb.pdb', 50)
	elif protocol == 'flxbb':
		RD.BDR(filename, 200)
		RD.flxbb('remodel.pdb', 50, 100)
		RD.Refine('flxbb.pdb', 50)
	else:
		print('Error in command string')

def MCProtocol(protocol, filename):
	RD = MCRosettaDesign()
	if protocol == 'fixbb':
		RD.fixbb(filename, 1.0, 10, 100, 'fixbb')
	elif protocol == 'flxbb':
		RD.flxbb(filename, 1.0, 10, 100, 'flxbb')
	else:
		print('Error in command string')

if __name__ == '__main__':
#	Protocol(sys.argv[1], sys.argv[2])
	MCProtocol(sys.argv[1], sys.argv[2])
