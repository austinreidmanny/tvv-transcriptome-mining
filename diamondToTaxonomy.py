#!/usr/bin/env python3

# Load the necessary libraries
import sys
import os
import requests

# Check to make sure this script is invoked correctly
if len(sys.argv) < 2:
	print("Usage {} {}".format(sys.argv[0], 'diamond_results.txt'))
	exit(1)

# Setup the inputs and outputs
inf = sys.argv[1]
outf = str(".".join(inf.split(".")[:-1])) + str(".taxonomy.txt")

# Check the first line of the Diamond file to make sure it's formatted correctly
with open(inf, "r") as input:
	for line in input:
		line = line.strip()
		line_arr = line.split()
		if len(line_arr) != 3:
			print("Diamond file should have the following format:")
			print("'{}\t{}\t{}'.format(sequence_name, taxonID, evalue)")
			exit(2)
		break                               

# Read in the DIAMOND results, find all the unique taxonIDs, and
# translate those into full taxonomic lineages
with open(inf, "r") as input:
	all_taxids = []

	# For each contig, if it belongs to a new taxon, store that taxon_id

	for line in input:
		line = line.strip()
		line_arr = line.split()	
		taxid = str(line_arr[1])
		if not taxid:
			continue
		if taxid not in all_taxids:
			all_taxids.append(taxid)
	# Now that all of the unique taxonomy ids have been read,
	# they can be translated into full taxonomy lineages

	taxid_translator = {}
	for taxid in all_taxids:

		# Take the TaxonID and query it against the taxonomy website
		url = 'http://taxonomy.jgi-psf.org/sc/simple/id/' + taxid
		t = requests.get(url)
		taxonomy = str(t.text)
		
		# Save that lineage with the taxid
		taxid_translator[taxid] = taxonomy

# Read the input file one more time, matching up the taxon_id with
# the taxonomy lineage that we now have for that taxon_id,
# and write all of that to an output file

with open(inf, "r") as input, \
     open(outf, "w") as output:
	for line in input:
		line = line.strip()
	
		# For each sequence, split it into an array and name each item
		line_arr = line.split()
		contig_name = str(line_arr[0])
		taxid = str(line_arr[1])
		evalue = str(line_arr[2])
		if not taxid:
			continue
	
		# Split the taxonomy result by rank
		taxonomy = taxid_translator[taxid]
		tax = {}
		tax_array = taxonomy.split(";")
		for pair in tax_array:
			rank = pair.split(":")
			type = rank[0]
			try:
				name = rank[1]
			except:
				name = 'N/A'
			tax[type] = name

		try:
			superkingdom = tax['sk']
		except:
			superkingdom = 'N/A'
		try:
			kingdom = tax['k']
		except:
			kingdom = 'N/A'
		try:
			phylum = tax['p']
		except:
			phylum = 'N/A'
		try:
			classs = tax['c']
		except:
			classs = 'N/A'
		try:
			order = tax['o']
		except:
			order = 'N/A'
	
		try:
			family = tax['f']
		except:
			family = 'N/A'
		try:
			genus = tax['g']
		except:
			genus = 'N/A'
		try:
			genus_species = tax['s']
			genus_species = genus_species.replace(" ", "_")
		except:
			genus_species = 'N/A'


		# Write it all to a new file
		output.write("\t".join((contig_name, taxid, evalue, \
                        	        superkingdom, kingdom, phylum, classs, \
                	                order, family, genus, genus_species)))
		output.write("\n")
