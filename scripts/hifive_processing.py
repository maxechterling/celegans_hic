#!/usr/bin/env python2.7

"""
Example:
	$ python hifive_processing.py alignments.raw name
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