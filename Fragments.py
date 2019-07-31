#!/usr/env/ python3

import os
import re
import sys
import bs4
import time
import datetime
import requests
import argparse
import urllib.request
from pyrosetta import *
from pyrosetta.toolbox import *
init('-out:level 0')

def Fragments(filename, username):
	pose = pose_from_pdb(filename)
	sequence = pose.sequence()
	web = requests.get('http://www.robetta.org/fragmentsubmit.jsp')
	payload = {	'UserName':				username,
				'Email':				'',
				'Notes':				'{}'.format(filename.split('.')[0]),
				'Sequence':				sequence,
				'Fasta':				'',
				'Code':					'',
				'ChemicalShifts':		'',
				'NoeConstraints':		'',
				'DipolarConstraints':	'',
				'type':					'submit'}
	session = requests.session()
	response = session.post('http://www.robetta.org/fragmentsubmit.jsp', data=payload, files=dict(foo='bar'))
	for line in response:
		line = line.decode()
		if re.search('<a href="(fragmentqueue.jsp\?id=[0-9].*)">', line):
			JobID = re.findall('<a href="(fragmentqueue.jsp\?id=[0-9].*)">', line)
	JobURL = 'http://www.robetta.org/' + JobID[0]
	ID = JobID[0].split('=')
	print('Fragments submitted to Robetta server [http://robetta.org/fragmentqueue.jsp]')
	print('Job ID: \u001b[32m{}\u001b[0m'.format(str(ID[1])))
	while True:
		Job = urllib.request.urlopen(JobURL)
		jobdata = bs4.BeautifulSoup(Job, 'lxml')
		status = jobdata.find('td', string='Status: ').find_next().text
		if status == 'Complete':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'), 'Status:', '\u001b[32m{}\u001b[0m'.format(status))
			break
		elif status == 'Active':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'), 'Status:', '\u001b[33m{}\u001b[0m'.format(status))
			time.sleep(180)
		else:
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M'), 'Status:', '\u001b[31m{}\u001b[0m'.format(status))
			time.sleep(300)
			continue
	sequence = pose.sequence()
	fasta = open('structure.fasta', 'w')
	fasta.write(sequence)
	fasta.close()
	time.sleep(1)
	os.system('wget http://www.robetta.org/downloads/fragments/'+str(ID[1])+'/aat000_03_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/'+str(ID[1])+'/aat000_09_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/'+str(ID[1])+'/t000_.psipred_ss2')
	os.rename('aat000_03_05.200_v1_3', 'frags.200.3mers')
	os.rename('aat000_09_05.200_v1_3', 'frags.200.9mers')
	os.rename('t000_.psipred_ss2', 'pre.psipred.ss2')
	frag = open('frags.200.9mers', 'r')
	data = open('RMSDvsPosition.dat', 'w')
	AVG = []
	for line in frag:
		if line.lstrip().startswith('position:'):
			line = line.split()
			size = line[1]
	frag.close()
	for i in range(1, int(size)):
		rmsd = []
		pose_copy = pyrosetta.Pose()
		pose_copy.assign(pose)
		frames = pyrosetta.rosetta.core.fragment.FrameList()
		fragset = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
		fragset.read_fragment_file('frags.200.9mers')
		fragset.frames(i, frames)
		movemap = MoveMap()
		movemap.set_bb(True)
		for frame in frames:
			for frag_num in range(1, frame.nr_frags()+1):
				frame.apply(movemap, frag_num, pose_copy)
				RMSD = rosetta.core.scoring.CA_rmsd(pose, pose_copy)
				rmsd.append(RMSD)
				lowest = min(rmsd)
				pose_copy.assign(pose)
		AVG.append(lowest)
		data.write(str(i)+'\t'+str(lowest)+'\n')
		print('\u001b[31mPosition:\u001b[0m {}\t\u001b[31mLowest RMSD:\u001b[0m {}\t|{}'.format(i, round(lowest, 3), '-'*int(lowest)))
	data.close()
	Average_RMSD = sum(AVG) / len(AVG)
	gnuplot = open('gnuplot_sets', 'w')
	gnuplot.write("""
	reset\n
	set terminal postscript\n
	set output './plot_frag.pdf'\n
	set encoding iso_8859_1\n
	set term post eps enh color\n
	set xlabel 'Position'\n
	set ylabel 'RMSD (\\305)'\n
	set yrange [0:]\n
	set xrange [0:]\n
	set xtics auto\n
	set xtics rotate\n
	set grid front\n
	unset grid\n
	set title 'Fragment Quality'\n
	set key off\n
	set boxwidth 0.5\n
	set style fill solid\n
	set label 'Average RMSD = {}' at graph 0.01 , graph 0.95 tc lt 7 font 'curior 12'\n
	plot 'RMSDvsPosition.dat' with boxes\n
	exit
	""".format(str(round(Average_RMSD, 3))))
	gnuplot.close()
	os.system('gnuplot < gnuplot_sets')
	os.remove('gnuplot_sets')
	print('\u001b[34mAverage RMSD:\u001b[0m {}'.format(round(Average_RMSD, 3)))
	return(Average_RMSD)

RMSD = Fragments(sys.argv[1], sys.argv[2])
