"""
python reg_mrna_ixns.py upreg.csv downreg.csv heatmap_dir
"""

import sys
import os
from math import log

import pandas as pandas
import numpy as np
import matplotlib.pyplot as plt

def rnaRegInput(file):
	chromDic = {'CHROMOSOME_I':'chrI', 'CHROMOSOME_II':'chrII', 'CHROMOSOME_III':'chrIII',\
				'CHROMOSOME_IV':'chrIV', 'CHROMOSOME_V':'chrV', 'CHROMOSOME_X':'chrX'}
	regDic = {}
	for value in chromDic.values():
		regDic[value] = []
	with open(file) as f:
		for line in f:
			if 'CHROMOSOME' in line:
				line = line.split(',')[1]
				chrom, coord = chromDic[line.split(':')[0]], line.split(':')[1]
				coord = (int(coord.split('-')[0]) + int(coord.split('-')[1])) / 2
				regDic[chrom].append(int(round(coord / 10000.) * 10000.) - 5000)
	return regDic

def getMeta():
	files = os.listdir(sys.argv[3])
	treatments, heatmaps = [], []
	for file in files:
		treatments.append(file.split('_')[0] + '_' + file.split('_')[1])
		if file.endswith('enrichment'):
			heatmaps.append(file)
	return list(set(heatmaps)), list(set(treatments))

def readMatrix(file):
	with open(sys.argv[3] + '/' + file) as f:
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
	chromDic = {'1':'chrI','2':'chrII','3':'chrIII','4':'chrIV', '5':'chrV', 'x':'chrX'}
	return midbins, np.array(matrix), chromDic[file.rstrip('.enrichment')[-1]]

def makeDistanceMatrix(midbins):
	matrix = []
	for i in range(len(midbins)):
		row = []
		for j in range(len(midbins)):
			row.append(abs(midbins[i] - midbins[j]))
		matrix.append(row)
	return np.array(matrix)

def scoreRegRegions(upregDic, downregDic, distMatrix, treatmentMatrices, fullbins):
	mindistance, maxdistance = 750000, 2500000
	upregScores, downregScores = [], []
	for chrom, matrix in treatmentMatrices.items():
		upreg, downreg = upregDic[chrom], downregDic[chrom]
		if len(matrix) != 0:
			for i in range(matrix[0].shape[1]):
				for j in range(matrix[0].shape[1]):
					if distMatrix[i][j] >= mindistance and distMatrix[i][j] <= maxdistance:
						if fullbins[i] in upreg and fullbins[j] in upreg:
							upregScores.append(matrix[0][i][j])
						elif fullbins[i] in downreg and fullbins[j] in downreg:
							downregScores.append(matrix[0][i][j])
	return [np.mean(upregScores), np.mean(downregScores)]

def plot(scoreDic):
	plt.figure()
	plt.bar(1, scoreDic['f1_wt'][0], color='r', label='F1 WT')
	plt.bar(2, scoreDic['f1_wt'][1], color='r')
	plt.bar(4, scoreDic['f1_morc1'][0], color='dimgray', label='F1 morc1(-/-)')
	plt.bar(5, scoreDic['f1_morc1'][1], color='dimgray')
	plt.bar(8, scoreDic['f4_wt'][0], color='tomato', label='F4 WT')
	plt.bar(9, scoreDic['f4_wt'][1], color='tomato')
	plt.bar(11, scoreDic['f4_morc1'][0], color='lightgray', label='F4 morc1(-/-)')
	plt.bar(12, scoreDic['f4_morc1'][1], color='lightgray')
	plt.show()
	plt.close()

def main():
	chroms = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrX']
	upregDic = rnaRegInput(sys.argv[1])
	downregDic = rnaRegInput(sys.argv[2])
	files, treatments = getMeta()
	fullbins, matrixDic = [], {}
	for treatment in treatments:
		matrixDic[treatment] = {}
		for chrom in chroms:
			matrixDic[treatment][chrom] = []
	for file in files:
		midbins, matrix, chrom = readMatrix(file)
		matrixDic[file.split('_')[0] + '_' + file.split('_')[1]][chrom].append(matrix)
		fullbins = fullbins + midbins
	distMatrix = makeDistanceMatrix(sorted(list(set(fullbins))))
	scoreDic = {}
	for treatment, treatmentMatrices in matrixDic.items():
		scoreDic[treatment] = scoreRegRegions(upregDic, downregDic, distMatrix, treatmentMatrices, fullbins)
	plot(scoreDic)

if __name__ == '__main__':
	main()