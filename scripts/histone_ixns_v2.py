"""
Example:
	$ python histone_ixns.py c_elegans_histone_mods.tsv heatmap_directory

This version plots genomic distance between loci versus contact enrichment between different
modification states.
"""

import sys
import os
from math import log
from math import isnan
from operator import itemgetter

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

def scoreTreatment(treatmentMatrices, distMatrix, stateDic):
	H3K36, H3K9, neither, both = {}, {}, {}, {}
	for dist in set(np.ndarray.flatten(distMatrix)):
		H3K36[dist], H3K9[dist], neither[dist], both[dist] = [], [], [], []
	for chrom, matrix in treatmentMatrices.items():
		if len(matrix) == 0:
			pass
		else:
			for i in range(matrix[0].shape[0]):
				for j in range(matrix[0].shape[0]):
					if stateDic[chrom][i] == 'H3K9me3' and stateDic[chrom][j] =='H3K9me3':
						H3K9[distMatrix[i][j]].append(matrix[0][i][j])
					elif stateDic[chrom][i] == 'H3K36me3' and stateDic[chrom][j] =='H3K36me3':
						H3K36[distMatrix[i][j]].append(matrix[0][i][j])
					elif stateDic[chrom][i] == 'neither' and stateDic[chrom][j] == 'neither':
						neither[distMatrix[i][j]].append(matrix[0][i][j])
	convertDic(H3K9)
	return convertDic(H3K9), convertDic(H3K36), convertDic(neither)

def convertDic(dic):
	outlist = []
	for dist, enrichments in dic.items():
		if isnan(np.mean(enrichments)) is False: 
			outlist.append([dist, np.mean(enrichments)])
		else:
			outlist.append([dist, 0.])
	#print outlist 
	outlist = sorted(outlist, key=itemgetter(0))
	x, y = [[x[0] for x in outlist], [x[1] for x in outlist]]
	lowresX, lowresY = [], []
	for i in range(0,len(x)-50,50):
		lowresX.append(np.mean(x[i:i+50]))
		lowresY.append(np.mean(y[i:i+50]))
	return lowresX, lowresY


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
	scoreDic = {}
	plt.figure()
	for treatment, treatmentMatrices in matrixDic.items():
		H3K9, H3K36, neither =scoreTreatment(treatmentMatrices, distMatrix, stateDic)
		plt.plot(H3K9[0][:-1], H3K9[1][:-1], label='%s H3K9' % (treatment))
		plt.plot(H3K36[0][:-1], H3K36[1][:-1], label='%s H3K36' % (treatment))
		plt.plot(neither[0][:-1], neither[1][:-1], label='%s neither' % (treatment))
	plt.legend()
	plt.show()
	plt.close()

if __name__ == '__main__':
	main()

