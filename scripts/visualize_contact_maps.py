"""
Example:
	$ python visualize_contact_maps.py hic.hdf5
"""

import sys

import hifive
import h5py

hic = hifive.HiC(sys.argv[1])

heatmap = hic.cis_heatmap(chrom='x',
					      binsize=50000,
						  arraytype='upper',
						  datatype='enrichment')

img = hifive.plotting.plot_upper_array(heatmap, symmetric_scaling=True)
img.save('f4_morc1_enr_heatmap_chrX.png')
