"""
Example:
	$ python histone_ixns.py c_elegans_histone_mods.tsv heatmap_directory
"""

import sys
import os
from math import log

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 	

def getMeta():
	files = os.listdir(sys.argv[2])
	treatments, heatmaps = [], []
	for file in files:
		treatments.append(file.split('_')[0] + '_' + file.split('_')[1])
		if file.endswith('enrichment'):
			heatmaps.append(file)
	return list(set(heatmaps)), list(set(treatments))

def readMatrix(file, tracksDF):
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
	stateTrack, chrom = findStates(midbins, file, tracksDF)
	return midbins, np.array(matrix), stateTrack, chrom

def findStates(midbins, file, tracksDF):
	track = []
	chrom = file.rstrip('.enrichment')[-1]
	chromDic = {'1':'chrI','2':'chrII','3':'chrIII','4':'chrIV', '5':'chrV', 'x':'chrX'}
	df = tracksDF.loc[tracksDF['chr'] == chromDic[chrom]]
	for midbin in midbins:
		state = df.loc[df['pos'] == midbin][['state']].values
		if len(state) == 0:
			track.append('NaN')
		else:
			track.append(state[0][0])
	return track, chromDic[chrom]

def makeDistanceMatrix(midbins):
	matrix = []
	for i in range(len(midbins)):
		row = []
		for j in range(len(midbins)):
			row.append(abs(midbins[i] - midbins[j]))
		matrix.append(row)
	return np.array(matrix)

def scoreTreatment(treatmentMatrices, distMatrix, stateDic, mindist, maxdist):
	H3K36, H3K9, neither, both = [], [], [], []
	for chrom, matrix in treatmentMatrices.items():
		if len(matrix) == 0:
			pass
		else:
			for i in range(matrix[0].shape[0]):
				for j in range(matrix[0].shape[0]):
					if distMatrix[i][j] >= mindist and distMatrix[i][j] <= maxdist:
						if stateDic[chrom][i] == 'H3K9me3' and stateDic[chrom][j] =='H3K9me3':
							H3K9.append(matrix[0][i][j])
						elif stateDic[chrom][i] == 'H3K36me3' and stateDic[chrom][j] =='H3K36me3':
							H3K36.append(matrix[0][i][j])
						elif stateDic[chrom][i] == 'neither' and stateDic[chrom][j] == 'neither':
							neither.append(matrix[0][i][j])
	return [np.mean(H3K9), np.std(H3K9), np.mean(H3K36), np.std(H3K36), np.mean(neither), np.std(neither)]

def plotting(scoreDic):
	plt.figure()
	for treatment, scores in scoreDic.items():
		if treatment == 'f1_wt':
			plt.bar(1, scores[0], yerr=scores[1], color='r', label='f1_wt')
			plt.bar(8, scores[2], yerr=scores[3], color='r')
			plt.bar(15, scores[4], yerr=scores[5], color='r')
		elif treatment == 'f4_wt':
			plt.bar(2, scores[0], yerr=scores[1], color='tomato', label='f4_wt')
			plt.bar(9, scores[2], yerr=scores[3], color='tomato')
			plt.bar(16, scores[4], yerr=scores[5], color='tomato')
		elif treatment == 'f1_morc1':
			plt.bar(4, scores[0], yerr=scores[1], color='dimgray', label='f1_morc1')
			plt.bar(11, scores[2], yerr=scores[3], color='dimgray')
			plt.bar(18, scores[4], yerr=scores[5], color='dimgray')
		elif treatment == 'f4_morc1':
			plt.bar(5, scores[0], yerr=scores[1], color='lightgray', label='f4_morc1')
			plt.bar(12, scores[2], yerr=scores[3], color='lightgray')
			plt.bar(19, scores[4], yerr=scores[5], color='lightgray')
	plt.xlabel('state of interacting fends')
	plt.ylabel('average interaction enrichment at short range (10kb-1Mb)')
	plt.legend()
	plt.savefig('histone_track_ixns.png')
	plt.close()


def main():
	chroms = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrX']
	tracksDF = pd.read_csv(sys.argv[1], sep='\t', header=None)
	tracksDF.columns = ['chr', 'pos', 'state']
	stateDic = {}
	files, treatments = getMeta()
	fullbins, matrixDic = [], {}
	for treatment in treatments:
		matrixDic[treatment] = {}
		for chrom in chroms:
			matrixDic[treatment][chrom] = []
	for file in files:
		midbins, matrix, stateTrack, chrom = readMatrix(file, tracksDF)
		stateDic[chrom] = stateTrack
		matrixDic[file.split('_')[0] + '_' + file.split('_')[1]][chrom].append(matrix)
		fullbins = fullbins + midbins
	distMatrix = makeDistanceMatrix(sorted(list(set(fullbins))))
	## Each treatment will have list with following values:
	## 0 H3K9 mean, 1 H3K9 stdev, 2 H3K36 mean, 3 H3K36 stdev, 4 neither mean, 5 neither stdev
	scoreDic = {}
	for treatment, treatmentMatrices in matrixDic.items():
		scoreDic[treatment] = scoreTreatment(treatmentMatrices, distMatrix, stateDic, 100000, 1000000)
	plotting(scoreDic)

if __name__ == '__main__':
	main()

#RK2GK-GNQWW-2RDQ4-JXHKR-R3KGJ