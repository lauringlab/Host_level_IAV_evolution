import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description=
	'This script takes a fasta file downloaded from GISAID and gathers \
	all the sequencing containing a label or one of a group of labels' ,usage ="")

parser.add_argument('fasta_file', metavar='fasta_file', nargs='+',
                    help='original fasta file')

parser.add_argument('output_file', metavar='output_file', nargs='+',
                    help='The output fasta file')

parser.add_argument('-OR', action='store',nargs="+",dest= 'OR',
                    help='The open reading frames we want to search for in the meta data of each sequence.\
                     separate multiple criteria with a space.')

                    


def main():
    args = parser.parse_args()

    fasta_file = args.fasta_file[0]
    output_file = args.output_file[0]
    identifier = args.OR
    if identifier==None:
    	raise ValueError, "Please provide an open reading frame"

    print "Reading sequences from %s writing to %s \n \
    looking for %s" %( fasta_file,output_file,identifier)

    # Get the desired Open Reading frames
    desired = []
    for record in SeqIO.parse(fasta_file,"fasta"):
    	if record.id in identifier:
    		desired.append(record)

    SeqIO.write(desired,output_file,"fasta")

if __name__ == '__main__':
    main()

