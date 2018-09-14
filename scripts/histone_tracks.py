"""
Example:
	python histone_tracks.py chip_seq_matrix chip_seq_series

"""

import sys

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def makeDF():
	df = pd.read_csv(sys.argv[1], sep='\t')
	with open(sys.argv[2]) as f:
		labelDic = {}
		for line in f:
			if line.startswith('!Sample_description'):
				for sample in line.split('\t'):
					if 'H3K36me3' in sample or 'H3K9me3' in sample:
						sample = sample.split()[2].split('_')
						labelDic[sample[-1].rstrip('"')] = '%s_%s_%s_%s' % \
							(sample[1], sample[0], sample[4], sample[2])
	labels = []
	for col in df.columns:
		if col in labelDic.keys():
			labels.append(labelDic[col])
		else:
			labels.append(col)
	df.columns = labels
	df = df[['chr', 'pos', 'F1_WT_H3K9me3_rep2', 'F1_WT_H3K36me3_rep2']]
	df.columns = ['chr', 'pos', 'f1_wt_H3K9me3', 'f1_wt_H3K36me3']
	return df

def plotThreshold(chrI_df):
	plt.figure()
	plt.plot(chrI_df[['pos']].head(n=1000), np.log(chrI_df[['f1_wt_H3K9me3']].head(n=1000)), label='H3K9me3')
	#plt.plot(chrI_df[['pos']].head(n=1000), np.log(chrI_df[['f1_wt_H3K36me3']].head(n=1000)), label='H3K36me3')
	plt.plot([-10000,1010000],[1.8,1.8], label='threshold @ 1.8')
	plt.xlabel('genomic position on chrI'), plt.ylabel('log(FPKM)')
	plt.title('H3K9me3 ChIP-Seq over chromosome I')
	plt.legend()
	plt.savefig('H3K9me3_track_chrI.png')
	#plt.show()

def assignState(df):
	states = []
	df.f1_wt_H3K36me3 = np.log(df['f1_wt_H3K36me3']) 
	df.f1_wt_H3K9me3 = np.log(df['f1_wt_H3K9me3'])
	for index, row in df.iterrows():
		if float(row['f1_wt_H3K36me3']) >= 2.3 and float(row['f1_wt_H3K9me3']) >= 1.8:
			states.append('both')
		elif float(row['f1_wt_H3K36me3']) >= 2.3:
			states.append('H3K36me3')
		elif float(row['f1_wt_H3K9me3']) >= 1.8:
			states.append('H3K9me3')
		else:
			states.append('neither')
	df['states'] = states
	return df[['chr','pos','states']]

def main():
	df = makeDF()
	## After plotting I've decided on the following thresholds on the LOG transformed data
	## H3K9me3 @ 1.8, H3K36me3 @ 2.3
	#plotThreshold(df.loc[df['chr'] == 'chrI'])
	stateDF = assignState(df)
	stateDF.to_csv('c_elegans_histone_mods.tsv', header=False, index=False, sep='\t')

if __name__ == '__main__':
	main()