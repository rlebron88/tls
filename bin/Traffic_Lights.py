#!/usr/bin/env python3

import os, sys, copy
import tabix
import numpy
import scipy.stats.mstats
#import statsmodels.sandbox.stats.multicomp

args = sys.argv[1].strip().split('\t')
chrom, start, end, strand, geneId, geneType, geneName = args[0:7]
start, end = int(start), int(end)
gene_exp = [float(tpm) if tpm != 'na' else 'na' for tpm in args[7:]]

out = open(os.path.join('Results', 'out', geneId + '.out'), 'wt')
err = open(os.path.join('Results', 'err', geneId + '.err'), 'wt')

cpgs = tabix.open('MethTable.sorted.gz')
upstream = 10**6
downstream = 10**6
#min_samples = 28

def meth_exp(meth, exp):
	while 'na' in meth:
		n = meth.index('na')
		meth.pop(n)
		exp.pop(n)
	while 'na' in exp:
		n = exp.index('na')
		meth.pop(n)
		exp.pop(n)
	return [meth, exp]

def cpgs_query(chrom, start, end, label):
	try:
		selected_cpgs = list(cpgs.query(chrom, start, end))
		if len(selected_cpgs) > 0:
			for cpg in selected_cpgs:
				cpg_chrom, cpg_start, cpg_end = cpg[0:3]
				meth = [float(mr) if mr != 'na' else 'na' for mr in cpg[3:]]
				global gene_exp
				exp = copy.deepcopy(gene_exp)
				meth, exp = meth_exp(meth, exp)
				if all([len(meth) == len(exp)]):
					corr, p = scipy.stats.mstats.spearmanr(meth, exp)
					line = [cpg_chrom, cpg_start, cpg_end, strand, geneId, geneType, geneName, corr, p, len(meth), label]
					line = '\t'.join(list(map(str, line))) + '\n'
					if numpy.isnan(corr):
						err.write(line)
					else:
						out.write(line)
				else:
					line = [cpg_chrom, cpg_start, cpg_end, strand, geneId, geneType, geneName, len(meth), len(exp), label, 'len_meth!=len_exp']
					line = '\t'.join(list(map(str, line))) + '\n'
					err.write(line)
	except:
		line = [chrom, start, end, strand, geneId, geneType, geneName, label, 'no_cpgs']
		line = '\t'.join(list(map(str, line))) + '\n'

cpgs_query(chrom, start, end, 'gene_body')
if strand == '+':
	up_start = start - upstream
	if up_start < 0:
		up_start = 0
	up_end = start - 1
	cpgs_query(chrom, up_start, up_end, 'upstream')
	down_start = end + 1
	down_end = end + downstream
	cpgs_query(chrom, down_start, down_end, 'downstream')
elif strand == '-':
        up_start = start - upstream
        if up_start < 0:
                up_start = 0
        up_end = start - 1
        cpgs_query(chrom, up_start, up_end, 'downstream')
        down_start = end + 1
        down_end = end + downstream
        cpgs_query(chrom, down_start, down_end, 'upstream')
else:
        up_start = start - upstream
        if up_start < 0:
                up_start = 0
        up_end = start - 1
        cpgs_query(chrom, up_start, up_end, 'nearby')
        down_start = end + 1
        down_end = end + downstream
        cpgs_query(chrom, down_start, down_end, 'nearby')

