"""
This code computes the average enrichment of the 

Usage:
	$ python mRNA_target_ixns.py targets.csv heatmap_dir
"""

import sys
import os

import pandas as pd

def getMeta():
	files = os.listdir(sys.argv[2])
	treatments, heatmaps = [], []
	for file in files:
		treatments.append(file.split('_')[0] + '_' + file.split('_')[1])
		if file.endswith('enrichment'):
			heatmaps.append(file)
	return list(set(heatmaps)), list(set(treatments))

def getTargetRNA(file):
	coordinates = []
	with open(file) as f:
		for i, line in enumerate(f):
			if i == 0:
				continue
			line = line.rstrip('\r\n').split(',')
			if line[0] == '':
				continue
			chrom = 'chr' + line[1].split(':')[0].split('_')[1]
			mid = int((int(line[1].split(':')[1].split('-')[0]) +
				  int(line[1].split(':')[1].split('-')[1]))	/ 2)
			coordinates.append([chrom, mid])
	return pd.DataFrame(coordinates, columns=['chr','coord'])

def main():
	heatmaps, treatments = getMeta()
	RNAcoords = getTargetRNA(sys.argv[1])
	print(RNAcoords)



if __name__ == '__main__':
	main()