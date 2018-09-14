#!/usr/bin/env python

"""
Python library for working with slurm
"""

import sys
import os 

class slurmScript(object):
	"""
	Creates a template slurm script with default arguments upon initialization. 
	User can specify variables for memory allocation, output/error files, etc
	and submit slurm jobs.

	Example:
		$ testScr = slurmScript()
		$ testScr.cmd = 'bowtie2 -x ./index -U reads.fastq -S reads.sam'
		$ testScr.modules = []
		$ testScr.writeScript()
	"""

	def __init__(self, cmd='No cmd', jobName='undefined_job', time='4:0:0',\
		         nodes=1, ntasks=1, partition='shared', rootpath='./', modules=[]):
		super(slurmScript, self).__init__()
		loadoutOptions = {'na':[], 'alignment':['bowtie2', 'misc']}
		self.cmd = cmd
		self.jobName = jobName
		## d-hr:min:sec
		self.time = time
		self.nodes = nodes
		self.ntasks = ntasks
		## e.g. shared, debug
		self.partition = partition
		self.rootpath = rootpath
		self.modules = modules
		self.writeScript()

	def __repr__(self):
		with open('%s%s.sh' % (self.rootpath, self.jobName)) as f:
			print(f.read())

	def writeScript(self):
		"""writes all specified variables to an output slurm script"""
		os.system('touch %s%s.sh' % (self.rootpath, self.jobName))
		with open('%s%s.sh' % (self.rootpath, self.jobName), 'w') as f:
			f.write('#!/bin/bash -l\n\n')
			f.write('#SBATCH\n')
			f.write('#SBATCH --job-name=%s\n' % (self.jobName))
			f.write('#SBATCH --time=%s\n' % (self.time))
			f.write('#SBATCH --error=%s%s\n' % (self.rootpath, self.jobName + '.err'))
			f.write('#SBATCH --output=%s%s\n' % (self.rootpath, self.jobName + '.out'))
			f.write('#SBATCH --nodes=%s\n' % (self.nodes))
			if len(self.modules) == 0:
				f.write('#SBATCH --ntasks-per-node=%s\n\n' % (self.ntasks))
			else:
				f.write('#SBATCH --ntasks-per-node=%s\n' % (self.ntasks))
				for mod in self.modules:
					f.write('module load %s' % (mod))
			f.write('\n' + self.cmd)

