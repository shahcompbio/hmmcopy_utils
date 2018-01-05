'''
Created on Aug 30, 2017

@author: dgrewal
'''


from fasta_to_read import FastaToRead
from bowtie_index import bowtieIndex
from map_counter import CountMappability
from read_to_map import readToMap

import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('reference',
                        help="path to the reference fasta file. The reference"\
                       " must be indexed with bowtie-build and samtools index")

    parser.add_argument('temp_dir',
                        help="path to a temporary dir for storing intermediate"\
                        " files. the directory will be removed after the run is"\
                        " complete if the clean up flag is set.")

    parser.add_argument('output',
                        help="path to the output mappability wig file")

    parser.add_argument('--chromosomes',
                        nargs='*',
                        default=map(str, range(1,23)) + ['X','Y'],
                        help="specify target chromosomes")

    parser.add_argument('--window_size',
                        default=50,
                        type=int,
                        help="specify window size for simulating reads")

    parser.add_argument('--mapcounter_window_size',
                        default=1000,
                        type=int,
                        help='specify the window size for mapcounter wig file'\
                        ' generation')


    parser.add_argument('--aligner',
                        default='bowtie',
                        choices=['bowtie'],
                        help='specify the aligner to use for aligning the'\
                        ' simulated reads')

    parser.add_argument('--maxhits',
                        default=4,
                        help='threshold for the num_hits value from bowtie')

    parser.add_argument('--cleanup',
                        default=False,
                        action='store_true',
                        help='if set, the tempdir will be deleted at the end'\
                        ' of execution')

    
    args = parser.parse_args()
    return args


class MapCounter(object):
    """
    runs all subprograms for mapcounter in one command
    keeps format consistent between read counter, mapcounter and read counter
    """
    
    def __init__(self, temp_dir, output, reference, chromosomes, aligner, maxhits, window_size,mapcounter_window_size, cleanup):
        self.temp_dir = temp_dir
        self.output = output
        self.reference = reference
        self.chromosomes = chromosomes
        self.aligner = aligner
        self.maxhits = maxhits
        self.window_size = window_size
        self.mapcounter_window_size = mapcounter_window_size
        self.cleanup = cleanup
    
    
    def main(self):
        fasta_read = os.path.join(self.temp_dir, 'fasta_reads.txt.gz')
        aln_out = os.path.join(self.temp_dir, 'aln_out.sam')
        bigwig = os.path.join(self.temp_dir, 'bigwig.txt')
    
    
        FastaToRead(self.reference, fasta_read, self.chromosomes, self.window_size).main()
        bowtieIndex(fasta_read, aln_out, self.reference, self.aligner).main()
        readToMap(aln_out, bigwig, self.reference, self.chromosomes, self.maxhits).main()
        CountMappability(bigwig, self.output, self.mapcounter_window_size).main()

        if self.cleanup:
            try:
                os.rmdir(self.temp_dir)
            except OSError:
                raise


if __name__ == '__main__':
    
    args = parse_args() 


    mapc = MapCounter(args.temp_dir, args.output, args.reference,
                      args.chromosomes, args.aligner, args.maxhits,
                      args.window_size, args.mapcounter_window_size,
                      args.cleanup)

    mapc.main()
