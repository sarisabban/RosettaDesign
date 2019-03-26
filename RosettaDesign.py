#!/usr/bin/python3

import os
import sys
import json
import Bio.PDB
import argparse
from pyrosetta import *
from pyrosetta.toolbox import *
init()

parser = argparse.ArgumentParser(description='RosettaDesign')
parser.add_argument('-i', '--fixbb', nargs='+', help='Performs fixed-backbone design')
parser.add_argument('-l', '--flxbb', nargs='+', help='Performs flexible-backbone design')
args = parser.parse_args()

class RosettaDesign(object):
	def __init__(self, filename):
		''' Generate the resfile '''
		self.filename = filename
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
		sasalist = []
		for x in dssp:
			if x[1] == 'A':
				sasa = 129*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'V':
				sasa = 174*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'I':
				sasa = 197*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'L':
				sasa = 201*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'M':
				sasa = 224*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'P':
				sasa = 159*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'Y':
				sasa = 263*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'F':
				sasa = 240*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'W':
				sasa = 285*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'R':
				sasa = 274*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'N':
				sasa = 195*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'C':
				sasa = 167*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'Q':
				sasa = 225*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'E':
				sasa = 223*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'G':
				sasa = 104*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'H':
				sasa = 224*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'K':
				sasa = 236*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'S':
				sasa = 155*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'T':
				sasa = 172*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'D':
				sasa = 193*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			if x[2] == 'G' or x[2] == 'H' or x[2] == 'I':	ss = 'H'
			elif x[2] == 'B' or x[2] == 'E':				ss = 'S'
			elif x[2] == 'S' or x[2] == 'T' or x[2] == '-':	ss = 'L'
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

	def __del__(self):
		''' Remove the resfile '''
		os.remove('resfile')
		try:	os.remove('fixbb.fasc')
		except:	os.remove('flxbb.fasc')

	def choose(self):
		''' Choose the lowest scoring structure '''
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
		try:	os.system('rm fixbb_*')
		except:	os.system('rm fixbb_*')

	def fixbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while maintaining a fixed backbone.
		Generates the structure.pdb file
		'''
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		packtask = standard_packer_task(pose)
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose, packtask, 'resfile')
		fixbb = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask, 10)
		job = PyJobDistributor('fixbb', 100, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			fixbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		self.choose()

	def flxbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while allowing for a flexible backbone.
		Generates the structure.pdb file
		'''
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		resfile = rosetta.core.pack.task.operation.ReadResfile('resfile')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		flxbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		flxbb.set_task_factory(task)
		flxbb.set_movemap(movemap)
		flxbb.set_scorefxn(scorefxn)
		job = PyJobDistributor('flxbb', 100, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			flxbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		self.choose()

def main():
	if args.fixbb:
		RD = RosettaDesign(sys.argv[2])
		RD.fixbb()
	elif args.flxbb:
		RD = RosettaDesign(sys.argv[2])
		RD.flxbb()

if __name__ == '__main__': main()
