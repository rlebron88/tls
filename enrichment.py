#!/usr/bin/env python3

import argparse
import lib.gtf

gencode = lib.gtf.Gtf('GTF/gencode.v28.chr_patch_hapl_scaff.annotation.gtf')
transcript = gencode.read_feature('transcript')
Transcripts = dict()

for t in transcript:
	chrom, start, end, strand = t[0], t[3], t[4], t[6]
	gene_type, gene_id, transcript_id = t[-1]['gene_type'], t[-1]['gene_id'], t[-1]['transcript_id']
	Transcripts[transcript_id] = {
		'chrom'       : chrom,
		'start'       : int(start),
		'end'         : int(end),
		'strand'      : strand,
		'gene_type'   : gene_type,
		'gene_id'     : gene_id,
		'exons'       : [],
		'introns'     : [],
		'5UTR'        : [],
		'3UTR'        : [],
		'UTR'         : [],
		'start_codon' : [],
		'stop_codon'  : []
	}

exon = gencode.read_feature('exon')

for e in exon:
	start, end = int(e[3]), int(e[4])
	transcript_id = e[-1]['transcript_id']
	Transcripts[transcript_id]['exons'].append([start, end])

for t in Transcripts:
	exons = Transcripts[t]['exons']
	if len(exons) > 1:
		for n in range(1, len(exons)):
			exon_a = exons[n-1]
			exon_b = exons[n]
			start = exon_a[0]+1
			end = exon_b[1]-1
			Transcripts[t]['introns'].append([start, end])
	n_exons = len(Transcripts[t]['exons'])
	n_introns = len(Transcripts[t]['introns'])
	if not n_exons == n_introns+1:
		print('{},{}'.format(n_exons, n_introns))

utr = gencode.read_feature('UTR')

for u in utr:
	start, end = int(u[3]), int(u[4])
	transcript_id = u[-1]['transcript_id']
	Transcripts[transcript_id]['UTR'].append([start, end])

n_coding = 0
n_errors = 0

start_codon = gencode.read_feature('start_codon')
stop_codon  = gencode.read_feature('stop_codon')

n_start = 0
n_stop  = 0

for c in start_codon:
	start, end = int(c[3]), int(c[4])
	transcript_id = c[-1]['transcript_id']
	Transcripts[transcript_id]['start_codon'].append([start, end])
	n_start += 1	

for c in stop_codon:
	start, end = int(c[3]), int(c[4])
	transcript_id = c[-1]['transcript_id']
	Transcripts[transcript_id]['stop_codon'].append([start, end])
	n_stop += 1

#exons_file   = open('GTF/exons.bed', 'wt')
#introns_file = open('GTF/introns.bed', 'wt')

for t in Transcripts:
	n = t
	t = Transcripts[t]
	gene_type = t['gene_type']
	print(gene_type)
	if any([gene_type == 'protein_coding', '_gene' in gene_type]):
		exons_file       = open('GTF/protein_coding_exons.bed', 'at')
		introns_file     = open('GTF/protein_coding_introns.bed', 'at')
		transcripts_file = open('GTF/protein_coding_transcripts.bed', 'at')
	else:
		exons_file       = open('GTF/non_coding_exons.bed', 'at')
		introns_file     = open('GTF/non_coding_introns.bed', 'at')
		transcripts_file = open('GTF/non_coding_transcripts.bed', 'at')
	chrom, start, end, name, score, strand = t['chrom'], t['start'], t['end'], n, '1000', t['strand']
	if int(start) > int(end):
		start, end = end, start
	line = list(map(str, [chrom, start, end, name, score, strand]))
	line = '\t'.join(line) + '\n'
	transcripts_file.write(line)
	exons, introns = t['exons'], t['introns']
	n_e = 1 if strand == '+' else len(exons)
	n_i = 1 if strand == '+' else len(introns)
	for exon in exons:
		e_start, e_end = exon
		if int(e_start) > int(e_end):
			e_start, e_end = e_end, e_start
		line = list(map(str, [chrom, e_start, e_end, '{}#exon{}'.format(name, n_e), score, strand]))
		line = '\t'.join(line) + '\n'
		exons_file.write(line)
		n_e += 1 if strand == '+' else -1
	for intron in introns:
		i_start, i_end = intron
		if int(i_start) > int(i_end):
			i_start, i_end = i_end, i_start
		line = list(map(str, [chrom, i_start, i_end, '{}#intron{}'.format(name, n_i), score, strand]))
		line = '\t'.join(line) + '\n'
		introns_file.write(line)
		n_i += 1 if strand == '+' else -1

