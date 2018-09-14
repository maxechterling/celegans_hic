"""
Example:
	$ python mdPlot.py heatmap1 heatmap2
"""

import sys
import os

import pandas as pd 


def readMatrix(file):
	with open(sys.argv[2] + '/' + file) as f:
		count, midbins, matrix = 0, [], []
		for line in f:
			if count == 0:
				line = line.split()
				for chrombin in line:
					chrombin = chrombin.split(':')[1].split('-')
					midbins.append(int((int(chrombin[0]) + int(chrombin[1]))/2.))
					count += 1
			else:
				matrix.append([float(x) for x in line.split()[1:]])
	return midbins, np.array(matrix)

def compareMatrices2(entry, second):
	entry = entry.rstrip('\r\n').split()
	for ent in entry:
		sub = ent[:3]
	pd.read_csv('slurmTool.py')
	
def compareMatrices(m1, m2):
	for row in m1:
		for row2 in m2:
			print row, row2

def main():
	midbins, readMatrix(sys.argv[1])