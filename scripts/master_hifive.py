"""
"""

import hifive
import sys

rawAlign, name = sys.argv[1], sys.argv[2] # Name will be the prefix of output files

## Load in the restriction enzyme digested fend coordinates
fend = hifive.Fend('%s_fend.hdf5' % (name), mode='w')
fend.load_fends('../ce10nm2.bed', genome_name='ce10', re_name='DpnII', format='bed')
fend.save()

## Load in the read data
data = hifive.HiCData('%s_data.hdf5' % (name), mode='w')
data.load_data_from_bam('%s_fend.hdf5' % (name), rawAlign, maxinsert=500)
data.save()

## Create a HiC object
hic = hifive.HiC('%s_hic.hdf5' % (name), 'w')
hic.load_data('%s_data.hdf5' % (name))
hic.save()

## Find HiC distance parameters
hic.find_distance_parameters(numbins=90, minsize=200, maxsize=0)

## Filter options. Chosen from filtering parameter sweep.
hic.filter_fends(mininteractions=10, mindistance=100000)
hic.save()

## Binning algorithm.
hic.find_binning_fend_corrections(max_iterations=1000,
                                  mindistance=5000000,
                                  maxdistance=0,
                                  num_bins=[20],
                                  model=['len'],
                                  parameters=['even'],
                                  usereads='cis',
                                  learning_threshold=1.0)

## Probability algorithm normalization
hic.find_probability_fend_corrections(mindistance=5000000,
                                      learningstep=0.5,
                                      max_iterations=1000,
                                      minchange=0.0005,
                                      precorrect=True)

hic.save()