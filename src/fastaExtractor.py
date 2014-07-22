#!/usr/bin/env python
'''
script to extract fasta sequences from original file by 
an index file containing the desired headers
'''
#@author: Philipp Sehnert
#@contact: philipp.sehnert[a]gmail.com

# IMPORTS
import sys, os
from argparse import ArgumentParser
from Bio import SeqIO

    
def main(argv = None):

    # Setup cmd interface
    parser = ArgumentParser(description = '%s -- extract sequences from fasta file by header' % 
                            (os.path.basename(sys.argv[0])),
                            epilog = 'created by Philipp Sehnert')
    parser.add_argument('--header', dest = 'header', required = True,
                        help = 'file with extracted headers')
    parser.add_argument('--fasta', dest = 'fasta', required = True,
                        help = 'fasta input file')
    parser.add_argument('--output',  dest = 'output', required = True,
                        help = 'fasta output file')
    args = parser.parse_args()

    if __name__ == '__main__':
    	sys.stdout.write('Creating index ...\n')
    	# build up a dictionary with the input of the header file
    	wanted = set()
    	with open(args.header, 'r') as header:
    			count = 0
    		# add header line by line to a dictionary
			for line in header:
				count += 1
				wanted.add(line.strip())
			sys.stdout.write('Loaded %d header to index\n' % (count))

			# load fasta sequences to memory
			fasta_sequences = SeqIO.parse(open(args.fasta,'r'), 'fasta')
			count = 0

			sys.stdout.write('Extracting sequences ...\n')
			with open(args.output, 'w') as output:
				for seq in fasta_sequences:
					if seq.id in wanted:
						count += 1
						SeqIO.write([seq], output, "fasta") 
    	
    	sys.stdout.write("Extraction of %d sequences successfull\n" % (count))    
sys.exit(main())
