
"""
Script for estimating distance function, filtering, and normalizing HiC data with hifive. Chains two \
normalization approaches: the binning algorithm followed by the probability algorithm.

Usage:
	$ ./hifive_post_processing.py input_hic output_hic
"""

import hifive
import sys

hic = hifive.HiC(sys.argv[1])

hic.find_distance_parameters(numbins=90, minsize=200, maxsize=0)

## Filter options. I'll have to experiment with these parameters.
hic.filter_fends(mininteractions=10, mindistance=100000)
hic.save()

## Binning algorithm
hic.find_binning_fend_corrections(max_iterations=1000,
                                      mindistance=5000000,
                                      maxdistance=0,
                                      num_bins=[10, 10],
                                      model=['len', 'gc'],
                                      parameters=['even', 'fixed-const'],
                                      usereads='cis',
                                      learning_threshold=1.0)

## Probability algorithm normalization
hic.find_probability_fend_corrections(mindistance=5000000,
                                      learningstep=0.5,
                                      max_iterations=1000,
                                      minchange=0.0005,
                                      precorrect=True)

hic.save()
