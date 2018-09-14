#!/usr/bin/env python

"""
"""

import sys

def readSams():
	for read1, read2 in zip(open(sys.argv[1]), open(sys.argv[2])):
		if read1.startswith('@'):
			pass
		else:
			read1, read2 = read1.rstrip('\r\n').split('\t'), read2.rstrip('\r\n').split('\t')
			## check mapping quality
			if int(read1[4]) >= 30 and int(read2[4]) >= 30:
				print('%s\t%s\t+\t%s\t%s\t+' % (convChrom(read1[2]), read1[3], convChrom(read2[2]), read2[3]))

def convChrom(chrom):
	chromDic = {'X':'x', 'I':'1', 'II':'2', 'III':'3', 'IV':'4', 'V':'5', 'M':'m'}
	return chromDic[chrom]



def main():
	readSams()

if __name__ == '__main__':
	main()