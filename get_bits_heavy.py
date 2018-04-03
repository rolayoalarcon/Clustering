# coding=utf-8
#		PROGRAM INFO.
#       Name: "Get bits".
#       Author: Vicente Fajardo.
#       Get a matrix of the largest bit score between pairs of sequences given as a fasta file.
#       Parameters:
#       --inputFile A fasta file with the sequences of interest.
#       --outputFiles
#				File to save the result of the blast search.
#				File to save the matrix.

#		LIBRARIES.
from re import split
from xml.dom import minidom
import pandas as pd
import subprocess as sub
import time

##		MAIN PROGRAM.

#	Carry out a blast search.
from Bio.Blast.Applications import NcbiblastpCommandline
blastp_cmd = NcbiblastpCommandline(query = "family-1.A.17.faa", subject = "family-1.A.17.faa", evalue = 500, outfmt = 5, max_hsps = 1, use_sw_tback = True, out = "family-1.A.17_blastp.xml")
stdout, stderr = blastp_cmd()

#	Parse the blast results in order to come up with the desired information.
#		Constrain: We are suppousing that there is an alignment result between each pair of sequences.
blast_records = minidom.parse("family-1.A.17_blastp.xml")

#	Let us define a blast record as the set of results (alignments) for a specific query.
blast_records = blast_records.getElementsByTagName("Iteration")

bits_blastpairs = {}
for blast_record in blast_records:
	query = blast_record.getElementsByTagName("Iteration_query-def")[0].firstChild.data
	bits_vals = []
	indices = []
	#bits_blastpairs[str(query)] = [[]]
	hits = blast_record.getElementsByTagName("Hit")
	for hit in hits:
		subject = hit.getElementsByTagName("Hit_id")[0].firstChild.data
		indices.append(str(subject))
		hsp = hit.getElementsByTagName("Hit_hsps")[0].getElementsByTagName("Hsp")[0]
		bits = hsp.getElementsByTagName("Hsp_bit-score")[0].firstChild.data
		bits_vals.append(float(bits))
	bits_blastpairs[str(query)] = pd.Series(bits_vals, index = indices)

#	Finally, we transform it into a pandas dataframe object.
bits_blastpairs = pd.DataFrame(bits_blastpairs)

queries = []
for query in bits_blastpairs:
	queries.append(query)

#	Keep the highest value for each pair.
for query1 in bits_blastpairs:
	for query2 in queries:
		bit1 = bits_blastpairs[query1][query2]
		bit2 = bits_blastpairs[query2][query1]
		bits_blastpairs[query1][query2] = max(bit1, bit2)
		bits_blastpairs[query2][query1] = max(bit1, bit2)

#	Write the result as a matrix to output file (csv format).
bits_blastpairs.to_csv(path_or_buf = "blastp_qands_bitscore_dist_matrix_protfam1A17.txt", sep = "\t", na_rep = "NA")

# Wait for file to write
time.sleep(5)

# Run analysis
rscript = 'Rscript Cluster_Protein_Analysis.R'
p = sub.Popen(rscript.split())