
"""
I need to determine the appropriate mininteraction and mindistance parameters for the HiC data
filtering. Since HiC hasn't really been done in C elegans I'm going to estimate these values with
a parameter sweep and see what percentage of fends are filtered out. One important thing to keep
in mind is that in other organisms TADs are on the order of ~1 Mb, so we can use that to hone in
on a good range for the mindistance parameter sweep.

Usage:
	python filtering_parameter_sweep.py hic_prefix
"""

import os
import sys

import hifive

## set parameter sweep range and step size
# low_ixn, high_ixn, step_ixn = 1, 50, 1
# low_dist, high_dist, step_dist = 1000, 10000000, 1000

# low_ixn, high_ixn, step_ixn = 1, 50, 2
# low_dist, high_dist, step_dist = 1000, 10000000, 10000

## medium resolution
low_ixn, high_ixn, step_ixn = 1, 50, 4
low_dist, high_dist, step_dist = 1000, 10000000, 100000


def HicSweep(prefix):
	## Redirect stdout in python so I can access the output later in the script
	# sys.stdout = open('%s.txt' % (prefix), 'w')
	for ixn in range(low_ixn, high_ixn + step_ixn, step_ixn):
		for dist in range(low_dist, high_dist + step_dist, step_dist):
			os.system('mkdir tmp_parameter_sweep')
			os.system('cp ../hifive/unprocessed_hic/%s* ./tmp_parameter_sweep/' % (prefix))
			hic = hifive.HiC('tmp_parameter_sweep/%s_hic.hdf5' % (prefix))
			print 'mininteraction=%s, mindistance=%s' % (ixn, dist)
			hic.filter_fends(mininteractions=ixn, mindistance=dist)
			os.system('rm -r tmp_parameter_sweep')
		
def main():
	HicSweep(sys.argv[1])

if __name__ == '__main__':
	main()
