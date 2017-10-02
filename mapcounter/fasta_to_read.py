'''
Created on Aug 28, 2017

@author: dgrewal
'''
import pysam
import argparse
from collections import deque
from itertools import islice


class FastaToRead(object):
    def __init__(self, ref_file, output, chromosomes, window_size):
        self.reference = ref_file
        self.window_size = window_size

        self.__get_fasta_reader()

        if chromosomes:
            self.chromosomes = chromosomes
        else:
            self.chromosomes = self.__get_chr_names()
        
        self.chr_lengths = self.__get_chr_lengths()
        
        self.output = output
        
    def __get_fasta_reader(self):
        self.fasta = pysam.FastaFile(self.reference)
    
    def __get_sequence(self, chrom, start, end):
        return self.fasta.fetch(chrom, start, end)

    def __get_chr_lengths(self):
        lengths = self.fasta.lengths
        names = self.fasta.references
        
        return {name:length for name,length in zip(names, lengths)}
    
    def __get_chr_names(self):
        return self.fasta.references
    
    
    def sliding_window(self, chrom, start, end):

        winsize = self.window_size * 2


        data = deque(self.__get_sequence(chrom, start, start+winsize))
 
        for pos in xrange(start, end):
 
            if self.window_size >= len(data):
                data.extend(self.__get_sequence(chrom, pos + self.window_size , pos + self.window_size + winsize))

            res = list(islice(data, 0, self.window_size))
            #remove one item from front
            data.popleft()
            
            yield pos, res

        
    
    def main(self):
        with open(self.output, 'w') as outfile:
            for chromosome in self.chromosomes:
                if chromosome not in self.chr_lengths:
                    raise Exception("Invalid chromosome")
            
                end = self.chr_lengths[chromosome] - self.window_size + 1
                assert end > 0
    
                for pos, seq in self.sliding_window(chromosome, 0, end):
                    seq = "".join(seq)
                    outfile.write("".join([">", chromosome,":", str(pos)]) + "\n")
                    outfile.write(seq+"\n")
            

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('reference')
    parser.add_argument('output')
    parser.add_argument('--chromosomes',
                        nargs='*')
    parser.add_argument('--window_size',
                        type=int,
                        default=50)

    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    args = parse_args()

    FastaToRead(args.reference, args.output, args.chromosomes, args.window_size).main()
