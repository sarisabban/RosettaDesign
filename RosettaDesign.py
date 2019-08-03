#!/usr/bin/ python3

import os
import sys
import json
import glob
import Bio.PDB
import argparse
from pyrosetta import *
from pyrosetta.toolbox import *
init('''
	-out:level 0
	-no_his_his_pairE
	-extrachi_cutoff 1
	-multi_cool_annealer 10
	-ex1 -ex2
	-use_input_sc
	''')

parser = argparse.ArgumentParser(description='RosettaDesign')
parser.add_argument('-i', '--fixbb', nargs='+', help='Performs fixed-backbone design')
parser.add_argument('-l', '--flxbb', nargs='+', help='Performs flexible-backbone design')
args = parser.parse_args()

class RosettaDesign(object):
	def __init__(self, filename):
		''' Generate the resfile. '''
		AminoAcid = {	'A':129, 'P':159, 'N':195, 'H':224,
						'V':174, 'Y':263, 'C':167, 'K':236,
						'I':197, 'F':240, 'Q':225, 'S':155,
						'L':201, 'W':285, 'E':223, 'T':172,
						'M':224, 'R':274, 'G':104, 'D':193}
		self.filename = filename
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
		sasalist = []
		for aa in dssp:
			sasa = AminoAcid[aa[1]]*aa[3]
			if sasa <= 25:      sasa = 'C'
			elif 25 < sasa < 40:sasa = 'B'
			elif sasa >= 40:    sasa = 'S'
			if aa[2] == 'G' or aa[2] == 'H' or aa[2] == 'I':   ss = 'H'
			elif aa[2] == 'B' or aa[2] == 'E':                 ss = 'S'
			elif aa[2] == 'S' or aa[2] == 'T' or aa[2] == '-': ss = 'L'
			sasalist.append((aa[0], aa[1], ss, sasa))
		resfile = open('.resfile', 'a')
		#resfile.write('NATRO\nEX 1\nEX 2\nUSE_INPUT_SC\n')
		resfile.write('START\n')
		for n, r, a, s in sasalist:
			if s == 'S' and a == 'L':   line = '{} A PIKAA PGNQSTDERKH\n'.format(n)
			elif s == 'S' and a == 'H': line = '{} A PIKAA EHKQR\n'.format(n)
			elif s == 'S' and a == 'S': line = '{} A PIKAA DEGHKNPQRST\n'.format(n)
			elif s == 'B' and a == 'L': line = '{} A PIKAA ADEFGHIKLMNPQRSTVWY\n'.format(n)
			elif s == 'B' and a == 'H': line = '{} A PIKAA ADEHIKLMNQRSTVWY\n'.format(n)
			elif s == 'B' and a == 'S': line = '{} A PIKAA DEFHIKLMNQRSTVWY\n'.format(n)
			elif s == 'C' and a == 'L': line = '{} A PIKAA AFGILMPVWY\n'.format(n)
			elif s == 'C' and a == 'H': line = '{} A PIKAA AFILMVWY\n'.format(n)
			elif s == 'C' and a == 'S': line = '{} A PIKAA FILMVWY\n'.format(n)
			resfile.write(line)
		resfile.close()
		self.SASA = sasalist
		# aa_composition file
		with open('.comp', 'w')as comp:
			comp.write("""
				PENALTY_DEFINITION
				PROPERTIES AROMATIC
				NOT_PROPERTIES POLAR CHARGED
				FRACTION 0.1
				PENALTIES 100 0 100
				DELTA_START -1
				DELTA_END 1
				BEFORE_FUNCTION CONSTANT
				AFTER_FUNCTION CONSTANT
				END_PENALTY_DEFINITION
				""")
		# netcharge file
		with open('.charge', 'w')as comp:
			comp.write("""
				DESIRED_CHARGE 0
				PENALTIES_CHARGE_RANGE -1 1
				PENALTIES 10 0 10
				BEFORE_FUNCTION QUADRATIC
				AFTER_FUNCTION QUADRATIC
				""")
		self.pose = pose_from_pdb(self.filename)
		# pushback aa_composition
		comp = pyrosetta.rosetta.protocols.aa_composition.AddCompositionConstraintMover()
		comp.create_constraint_from_file('.comp')
		comp.apply(self.pose)
		# pushback netcharge
		charge = pyrosetta.rosetta.protocols.aa_composition.AddNetChargeConstraintMover()
		charge.create_constraint_from_file('.charge')
		charge.apply(self.pose)
		self.starting_pose = Pose()
		self.starting_pose.assign(self.pose)
		self.scorefxn = get_fa_scorefxn()
		self.scorefxn_G = get_fa_scorefxn()
		AAcomp      = pyrosetta.rosetta.core.scoring.ScoreType.aa_composition
		NETq        = pyrosetta.rosetta.core.scoring.ScoreType.netcharge
		AArep       = pyrosetta.rosetta.core.scoring.ScoreType.aa_repeat
		ASPpen      = pyrosetta.rosetta.core.scoring.ScoreType.aspartimide_penalty
		HBnet       = pyrosetta.rosetta.core.scoring.ScoreType.hbnet
		MHCep       = pyrosetta.rosetta.core.scoring.ScoreType.mhc_epitope
		VOIDpen     = pyrosetta.rosetta.core.scoring.ScoreType.voids_penalty
		ABurUnsatPen= pyrosetta.rosetta.core.scoring.ScoreType.approximate_buried_unsat_penalty
		BurUnsatPen = pyrosetta.rosetta.core.scoring.ScoreType.buried_unsatisfied_penalty
		#self.scorefxn_G.set_weight(AAcomp,      1.00)
		#self.scorefxn_G.set_weight(NETq,        1.00)
		#self.scorefxn_G.set_weight(HBnet,       1.00)
		#self.scorefxn_G.set_weight(VOIDpen,     0.10)
		self.scorefxn_G.set_weight(AArep,       1.00)
		self.scorefxn_G.set_weight(ASPpen,      1.00)
		self.scorefxn_G.set_weight(MHCep,       0.00)
		self.scorefxn_G.set_weight(BurUnsatPen, 1.00)
		self.scorefxn_G.set_weight(ABurUnsatPen,5.00)
		self.relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		self.relax.set_scorefxn(self.scorefxn)
	def __del__(self):
		''' Remove the resfile. '''
		os.remove('.resfile')
		os.remove('.comp')
		os.remove('.charge')
		for f in glob.glob('f[il]xbb.fasc'): os.remove(f)
	def choose(self):
		''' Choose the lowest scoring structure. '''
		try:	scorefile = open('fixbb.fasc', 'r')
		except:	scorefile = open('flxbb.fasc', 'r')
		score = 0
		name = None
		for line in scorefile:
			line = json.loads(line)
			score2 = line.get('total_score')
			if score2 < score:
				score = score2
				name = line.get('decoy')
		os.system('mv {} structure.pdb'.format(name))
		for f in glob.glob('f[il]xbb_*'): os.remove(f)
	def fixbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while maintaining a fixed backbone.
		Generates the structure.pdb file.
		'''
		resfile = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('.resfile')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(False)
		movemap.set_chi(True)
		fixbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		fixbb.set_task_factory(task)
		fixbb.set_movemap(movemap)
		fixbb.set_scorefxn(self.scorefxn_G)
		self.relax.apply(self.pose)
		job = PyJobDistributor('fixbb', 100, self.scorefxn)
		job.native_pose = self.starting_pose
		while not job.job_complete:
			self.pose.assign(self.starting_pose)
			fixbb.apply(self.pose)
			self.relax.apply(self.pose)
			job.output_decoy(self.pose)
	def flxbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while allowing for a flexible backbone.
		Generates the structure.pdb file.
		'''
		resfile = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('.resfile')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		flxbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		flxbb.set_task_factory(task)
		flxbb.set_movemap(movemap)
		flxbb.set_scorefxn(self.scorefxn_G)
		self.relax.apply(self.pose)
		job = PyJobDistributor('flxbb', 100, self.scorefxn)
		job.native_pose = self.starting_pose
		while not job.job_complete:
			self.pose.assign(self.starting_pose)
			flxbb.apply(self.pose)
			self.relax.apply(self.pose)
			job.output_decoy(self.pose)

def main():
	if args.fixbb:
		RD = RosettaDesign(sys.argv[2])
		RD.fixbb()
	elif args.flxbb:
		RD = RosettaDesign(sys.argv[2])
		RD.flxbb()

if __name__ == '__main__': main()
