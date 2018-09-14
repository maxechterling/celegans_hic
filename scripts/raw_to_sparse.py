#!/usr/bin/env python3

"""
Usage:
	./raw_to_sparse.py hic_data.raw chrom1,chrom2,...
"""

import sys

def makeChromDic():
	chromDic, chroms = {}, sys.argv[2].split(',')
	for chrom in chroms:
		chromDic[chrom] = {}
	print(chromDic)

	with open(sys.argv[1]) as f:
		for i, line in enumerate(f):
			print(line)
			line = line.rstrip('\r\n').split('\t')
			if line[0] in chromDic.keys():
				else:
					ixnKey = ','.join(ixnKey)
				if ixnKey in chromDic[line[0]].keys():
					chromDic[line[0]][ixnKey] += 1
				else: 
					chromDic[line[0]][ixnKey] = 1
				print(i)
	return chromDic

def writeOutput(chromDic):
	for chrom, ixnFreqs in chromDic.items():
		filename = sys.argv[1].split('.')[0] + '.chr' + chrom
		with open(filename, 'w') as f:
			f.write('region1\tregion2\tIF\n')
			for regions, freq in ixnFreqs.items():
				regions = regions.split(',')
				out = '%s\t%s\t%s\n' % (regions[0], regions[1], freq)
				f.write(out)

def main():
	chromDic = makeChromDic()
	writeOutput(chromDic)

if __name__ == '__main__':
	main()