#!/usr/bin/env python

"""
This version of samToRaw accounts for strand of the alignment
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
				if int(read1[1]) & 16:
					strand1 = '-'
				else:
					strand1 = '+'
				if int(read2[1]) & 16:
					strand2 = '-'
				else:
					strand2 = '+'
				print('%s\t%s\t%s\t%s\t%s\t%s' % (convChrom(read1[2]), read1[3], strand1,\
												  convChrom(read2[2]), read2[3], strand2))

def convChrom(chrom):
	chromDic = {'X':'x', 'I':'1', 'II':'2', 'III':'3', 'IV':'4', 'V':'5', 'M':'m'}
	return chromDic[chrom]



def main():
	readSams()

if __name__ == '__main__':
	main()