#!/usr/bin/env python3

import sys, requests, ujson

experiment = sys.argv[1]

r = requests.get('https://www.encodeproject.org/experiments/{}/?format=json'.format(experiment))
if r.status_code != 200:
	sys.stderr.write('{}: {}\n'.format(r.status_code, experiment))
else:
	data = r.json()['replicates']
	accession = set()
	donor = set()
	biosample = set()
	biosample_json = set()
	for replicate in data:
		accession.add(replicate['library']['biosample']['accession'])
		donor.add(replicate['library']['biosample']['donor']['uuid'])
		biosample.add(replicate['library']['biosample']['biosample_term_id'])
		biosample_json.add(ujson.dumps(replicate['library']['biosample']))
	accession = ';'.join(list(accession))
	donor = ';'.join(list(donor))
	biosample = ';'.join(list(biosample))
	biosample_json = ';'.join(list(biosample_json))
	output = '\t'.join([experiment, accession, donor, biosample, biosample_json]) + '\n'
	sys.stdout.write(output)

