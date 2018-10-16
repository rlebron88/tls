#!/usr/bin/env python3

import os.path, gzip, bz2, csv, collections

class Gtf:
	def __init__(self, f):
		ext = os.path.splitext(f)[1]
		if ext == '.gz':
			self.handle = gzip.open(f, 'rt')
		elif ext == '.bz2':
			self.handle = bz2.open(f, 'rt')
		else:
			self.handle = open(f, 'rt')
		self.reader = csv.reader(self.handle, delimiter = '\t')
	def read_feature(self, feature):
		self.handle.seek(0)
		for line in self.reader:
			if not line[0].startswith('#'):
				if line[2] == feature:
					additional_info = [field.strip().replace('"', '').split(' ', 1) for field in line[-1].strip().split(';')]
					line[-1] = collections.OrderedDict([tuple([field[0], int(field[1]) if field[1].isnumeric() else field[1]]) for field in additional_info if len(field) == 2])
					yield line

class GeneHancer:
	def __init__(self, f):
		ext = os.path.splitext(f)[1]
		if ext == '.gz':
			self.handle = gzip.open(f, 'rt')
		elif ext == '.bz2':
			self.handle = bz2.open(f, 'rt')
		else:
			self.handle = open(f, 'rt')
		self.header = self.handle.readline().strip().split('\t')
		self.reader = csv.reader(self.handle, delimiter = '\t')
	def read_feature(self, feature):
		self.handle.seek(0)
		self.handle.readline()
		for line in self.reader:
			if not line[0].startswith('#'):
				if line[2] == feature:
					additional_info = [field.strip().split('=') for field in line[-1].strip().split(';')]
					attributes = [additional_info[0], ['connected_genes', []]]
					for i in range(1, len(additional_info), 2):
						attributes[1][1].append(collections.OrderedDict({
							'gene_name' : additional_info[i][1],
							'score' : float(additional_info[i+1][1])
						}))
					line[-1] = collections.OrderedDict([[field[0], field[1]] for field in attributes if len(field) == 2])
					yield line
