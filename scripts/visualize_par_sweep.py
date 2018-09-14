"""
Usage:
	python visualize_par_sweep.py sweep.results sweep.par
"""

import sys

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def readSweep(sweepResults, sweepPars):
	filtResults, parameters = [], []
	with open(sweepResults) as f1:
		for line in f1:
			line = line.split(' ')
			print line
			totalfends, fxnRemain = int(line[5]), 1 - (float(line[3]) / float(line[5]))
			filtResults.append(fxnRemain)
	with open(sweepPars) as f2:
		for line in f2:
			line = line.rstrip('\r\n').split(',')
			print line
			ixn, dist = line[0].lstrip('mininteraction='), line[1].lstrip(' mindistance=')
			parameters.append((int(ixn), int(dist)))
	return parameters, filtResults, totalfends

def findDimensions(parameters):
	ixns, dists = [], []
	[ixns.append(each[0]) for each in parameters], [dists.append(each[1]) for each in parameters]
	return len(set(ixns)), len(set(dists))

def convertMatrix(parameters, filtResults):
	ixnLen, distLen = findDimensions(parameters)
	print ixnLen, distLen
	print len(filtResults)
	matrix = []
	for i in range(0, len(filtResults), distLen):
		matrix.append(filtResults[i:i+distLen])
	print matrix
	return matrix

def unique(sequence):
    seen = set()
    return [x for x in sequence if not (x in seen or seen.add(x))]

def getLabels(parameters):
	ixns, dists = [],[]
	for par in parameters:
		ixns.append(par[1]), dists.append(par[0])
	return unique(ixns), unique(dists)

def plotHeatmap(parameters, parameterMatrix, title):
	plt.figure()
	ylab, xlab = getLabels(parameters)
	ax = sns.heatmap(parameterMatrix)
	plt.xticks([i for i in range(0, len(ylab), 7)], ylab[0::7]), plt.yticks([i for i in range(0, len(xlab), 3)], xlab[0::3])
	plt.title(title)
	plt.xlabel('mindistance'), plt.ylabel('mininteractions')
	plt.tight_layout()
	plt.savefig(title.split(';')[0])

def main():
	parameters, filtResults, totalfends = readSweep(sys.argv[2], sys.argv[1])
	parameterMatrix = convertMatrix(parameters, filtResults)
	plotHeatmap(parameters, parameterMatrix, 'f4_hrde1; %s starting fends' % (totalfends))

if __name__ == '__main__':
	main()

# def readSweep(sweepFile):
# 	parameters, filtResults = [], []
# 	with open(sweepFile) as f:
# 		for line in f:
# 			if line.startswith('min') and line.endswith('fends\n'):
# 				continue
# 				#parameters.append('None'), filtResults.append(0)
# 			elif line.startswith('min'):
# 				line = line.rstrip('\r\n').split(',')
# 				ixn, dist = line[0].lstrip('mininteraction='),\
# 					        line[1].lstrip(' mindistance=')
# 				parameters.append((int(ixn), int(dist)))
# 			elif line.startswith('Filtering'):
# 				line = line.split(' ')
# 				totalfends, fxnRemain = int(line[5]), 1 - (float(line[3]) / float(line[5]))
# 				filtResults.append(fxnRemain)
# 	return parameters, filtResults, totalfends