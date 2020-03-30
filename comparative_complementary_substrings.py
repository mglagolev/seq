#!/usr/bin/python

from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
import argparse

parser = argparse.ArgumentParser(description = 'Find sequences which maintain complementary structures through mutation')

parser.add_argument('--files', type = str, nargs = '+', help = 'sequence files')

parser.add_argument('--min', type = int, nargs = '?', help = 'minimum substring length', default = 10)

parser.add_argument('--max', type = int, nargs = '?', help = 'maximum substring length', default = 20)

args = parser.parse_args()

#TODO:
#create 2d array
#upon finding coincidence update index (source, match)
#make 2d arrays for each of the sequences
#compare the matches in different sequences

def trim_seq(seq):
	seq = seq.replace("\n","")
	seq = seq.replace(" ","")
	seq = seq.replace("	","")
	return seq
	

#https://stackoverflow.com/questions/4664850/how-to-find-all-occurrences-of-a-substring
def find_all(a_str, sub, start = 0):
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += 1

def find_complementary_subseqs(seq, length, heat_matrix):
	all_indices = []
	s_seq = str(seq)
	s_comp_seq = str(seq.complement())
	for i in range(len(s_seq) - length):
		sub = s_seq[i:i+length]
		for index in find_all(s_comp_seq, sub, i+1)
			for k in range(length):
				for l in range(length):
					heat_matrix[i+k][index+l] += 1
	return
all_heat_seqs = {}

for length in range(args.min,args.max + 1):

	seqs = []

	for sfile in args.files:
		infile = open(sfile,'rb')
		seqs.append(Seq(trim_seq(infile.read()), generic_rna))
		infile.close()

	heat_range = max(map(len, seqs))
	heat_matrix = [[0 for x in range(heat_range)] for y in range(heat_range)]

	for seq in seqs:
		heat_matrix = [[0 for x in range(heat_range)] for y in range(heat_range)]
		find_complementary_subseqs(seq,length,heat_matrix)
		for indices in all_indices:
			for index in indices:
				for i in range(index, index+length+1):
					heat_seq[i] += 1
			
	all_heat_seqs[length] = heat_seq

for i in range(len(all_heat_seqs[args.min])):
	print i+1,
	for j in range(args.min, args.max + 1):
		print all_heat_seqs[j][i],
	print "\n",






