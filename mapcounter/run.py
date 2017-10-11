'''
Created on Aug 30, 2017

@author: dgrewal
'''


from fasta_to_read import FastaToRead
from bowtie_index import bowtieIndex
from mapcounter import MapCounter
from read_to_map import readToMap

import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('reference')
    parser.add_argument('temp_dir')
    parser.add_argument('output')
    parser.add_argument('--chromosomes',
                        nargs='*')

    parser.add_argument('--window_size',
                        default=50)

    parser.add_argument('--mapcounter_window_size',
                        default=1000)


    parser.add_argument('--aligner',
                        default='bowtie')

    parser.add_argument('--maxhits',
                        default=4)

    parser.add_argument('--cleanup',
                        default=False,
                        action='store_true')

    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    
    args = parse_args() 

    fasta_read = os.path.join(args.temp_dir, 'fasta_reads.txt')
    aln_out = os.path.join(args.temp_dir, 'aln_out.sam')
    bigwig = os.path.join(args.temp_dir, 'bigwig.txt')


    FastaToRead(args.reference, fasta_read, args.chromosomes, args.window_size).main()
    bowtieIndex(fasta_read, aln_out, args.reference, args.aligner).main()
    readToMap(aln_out, bigwig, args.maxhits).main()
    MapCounter(bigwig, args.output, args.mapcounter_window_size).main()

    if args.cleanup:
        try:
            os.rmdir(args.tempdir)
        except OSError:
            raise
