#!/usr/bin/env python3

import tabix
import scipy.stats.mstats
import statsmodels.sandbox.stats.multicomp
import sys

genes = open(sys.argv[1], 'rt')
cpgs = tabix.open(sys.argv[2])
tmp = open(sys.argv[3], 'wt')
output = open(sys.argv[4], 'wt')

pvalues = []

for line in genes:
	line = line.strip().split('\t')
	ID, chrom, start, end = line[0:4]
	data = line[4:]
	records = cpgs.query(chrom, int(start), int(end))
	for record in records:
		e = data.copy()
		chrom, start, end = record[0:3]
		d = record[3:]
		while 'N' in e:
			n = e.index('N')
			e.pop(n)
			d.pop(n)
		while 'N' in d:
			n = d.index('N')
			d.pop(n)
			e.pop(n)
		exp = [float(n) for n in e]
		meth = [float(n) for n in d]
		if len(exp) != len(meth):
			sys.stderr.write('exp: {}\nmeth: {}\n'.format(exp, meth))
		if len(exp) >= 17:
			scc = scipy.stats.mstats.spearmanr(exp, meth)
			corr = str(scc.correlation)
			p = str(float(scc.pvalue))
			pvalues.append(float(p))
			E, M = ','.join(map(str, exp)), ','.join(map(str, meth))
			new = '\t'.join([chrom, start, ID, E, M, corr, p]) + '\n'
			tmp.write(new)

genes.close()
tmp.close()

tmp = open(sys.argv[3], 'rt')

fdr = statsmodels.sandbox.stats.multicomp.multipletests(pvalues, alpha = 0.05, method = 'fdr_bh')
acepted = list(list(fdr)[0])
pvals_corrected = list(list(fdr)[1])

n = 0
for line in tmp:
	line = line.strip().split('\t')
	line = line + [str(acepted[n]), str(pvals_corrected[n])]
	line = '\t'.join(line) + '\n'
	output.write(line)
	n += 1

