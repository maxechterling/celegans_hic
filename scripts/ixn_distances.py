"""
Example:
	python ixn_distances.py wt_heatmap.txt morc_heatmap.txt
"""

import sys

import numpy as np
from operator import itemgetter
import matplotlib.pyplot as plt

def readMatrix(filename):
	with open(filename) as f:
		count, midbins, matrix = 0, [], []
		for line in f:
			if count == 0:
				line = line.split()
				for chrombin in line:
					chrombin = chrombin.split(':')[1].split('-')
					midbins.append(int((int(chrombin[0]) + int(chrombin[1]))/2.))
					count += 1
			else:
				matrix.append([int(x) for x in line.split()[1:]])
	return midbins, np.array(matrix)
				
def makeDistanceMatrix(midbins):
	matrix = []
	for i in range(len(midbins)):
		row = []
		for j in range(len(midbins)):
			row.append(abs(midbins[i] - midbins[j]))
		matrix.append(row)

	return np.array(matrix)

def contactProbability(observed, distMatrix):
	scoreDic = {}
	for dist in set(distMatrix.flatten()):
		scoreDic[dist] = 0
	for i in range(len(observed[0])):
		for j in range(len(observed[0])):
			scoreDic[distMatrix[i][j]] += observed[i][j]
	counts = sorted(scoreDic.iteritems(), key=itemgetter(0))[1:]
	return [np.log(each[0]) for each in counts], [np.log(float(each[1])/np.sum(observed)) for each in counts]

def plot(x, y, x1, y1, x2, y2, x3, y3):
	plt.figure()
	plt.plot(x, y, label='f1_wt')
	plt.plot(x1, y1, label='f1_morc1-/-')
	plt.plot(x2, y2, label='f4_wt')
	plt.plot(x3, y3, label='f4_morc1-/-')
	plt.legend()
	plt.xlabel('log distance(bp)')
	plt.ylabel('log contact probability')
	plt.title('Contact probability vs distance; chromosome x')
	plt.savefig('all_treatment_ixn_dist.png')
	plt.close()

def main():
	midbins, observed = readMatrix(sys.argv[1])
	distMatrix = makeDistanceMatrix(midbins)
	midbins2, observed2 = readMatrix(sys.argv[2])
	midbins3, observed3 = readMatrix(sys.argv[3])
	midbins4, observed4 = readMatrix(sys.argv[4])
	x, y = contactProbability(observed, distMatrix)
	x1, y1 = contactProbability(observed2, distMatrix)
	x2, y2 = contactProbability(observed3, distMatrix)
	x3, y3 = contactProbability(observed4, distMatrix)
	plot(x, y, x1, y1, x2, y2, x3, y3)

if __name__ == '__main__':
	main()