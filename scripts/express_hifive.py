import hifive
import sys

hic = hifive.HiC(sys.argv[1])

hic.find_distance_parameters(numbins=90, minsize=200, maxsize=0)

## Filter options. I'll have to experiment with these parameters.
hic.filter_fends(mininteractions=10, mindistance=100000)
hic.save()

hic.find_express_fend_corrections(iterations=1000,
								  mindistance=0,
								  usereads='cis',
								  remove_distance=True)

hic.save()