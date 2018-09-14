"""
Plot the average enrichment for all possible cis interactions at each distance across all chromosomes.
You need to create heatmaps with make_heatmap.py for HiC data before running this script. 

Heatmap directory must contain only the heatmaps for
chromosomes and samples you wish to include in the final plot.

Usage:
	$ python avg_enr_vs_dist.py heatmap_directory

"""

import sys
import os

import numpy as np 
import matplotlib.pyplot as plt 

def getMeta():
	files = os.listdir(sys.argv[1])
	treatments, enrichmentHeatmaps = [], []
	for file in files:
		treatments.append(file.split('_')[0] + '_' + file.split('_')[1])
		if file.endswith('enrichment'):
			enrichmentHeatmaps.append(file)
	return list(set(enrichmentHeatmaps)), list(set(treatments))

def readMatrix(file):
	with open(sys.argv[1] + file) as f:
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

def makeDistanceMatrix(midbins):
	matrix = []
	for i in range(len(midbins)):
		row = []
		for j in range(len(midbins)):
			row.append(abs(midbins[i] - midbins[j]))
		matrix.append(row)
	return np.array(matrix)

def scoreMatrices(matrices, distMatrix):
	enrichmentDic, coords = {}, []
	fullbins = list(set(np.ndarray.flatten(distMatrix)))
	for each in fullbins:
		enrichmentDic[each] = []
	for matrix in matrices:
		for i in range(matrix.shape[0]):
			for j in range(matrix.shape[0]):
				enrichmentDic[int(distMatrix[i][j])].append(matrix[i][j])
	for binsize, enrichments in enrichmentDic.items():
		coords.append((binsize, np.average(enrichments)))
	coords = sorted(coords, key=lambda x: x[0])
	return [np.log(x[0]) for x in coords], [x[1] for x in coords]

def avgEnrichment(distMatrix, matrixDic):
	plt.figure()
	for key, matrices in matrixDic.items():
		x, y = scoreMatrices(matrices, distMatrix)
		plt.plot(x[1:], y[1:], label=key)
	plt.xlabel('distance(bp)'), plt.ylabel('avg enrichment')
	plt.legend()
	plt.title('Average interaction enrichment vs distance')
	plt.show()
	plt.close()

def main():
	files, treatments = getMeta()
	fullbins, matrixDic = [], {}
	for treatment in treatments:
		matrixDic[treatment] = []
	for file in files:
		midbins, matrix = readMatrix(file)
		matrixDic[file.split('_')[0] + '_' + file.split('_')[1]].append(matrix)
		fullbins = fullbins + midbins
	distMatrix = makeDistanceMatrix(sorted(list(set(fullbins))))
	avgEnrichment(distMatrix, matrixDic)
	

if __name__ == '__main__':
	main()