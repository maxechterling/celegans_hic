"""
Usage:
	$ python make_heatmap.py hic.hdf5
"""

import sys
import hifive

hic = hifive.HiC(sys.argv[1])

hic.write_heatmap(filename='f4_morc1_10kb_enr', binsize=10000, datatype='enrichment', format='txt', includetrans=False)
hic.save()