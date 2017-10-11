'''
Created on Aug 28, 2017

@author: dgrewal
'''
import pysam
import argparse
from collections import deque
from itertools import islice


class FastaToRead(object):
    """simulate reads from the fasta file
    """

    def __init__(self, ref_file, output, chromosomes, window_size):
        """
        :param ref_file: reference fasta file
        :param output: output fastq file with simulated reads
        :param chromosomes: target chromosomes
        :param window_size: (int) window_size
        """
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
        """returns pysam fasta object
        :returns pysam fasta object
        """
        self.fasta = pysam.FastaFile(self.reference)
    
    def __get_sequence(self, chrom, start, end):
        """returns sequence from the specified region
        :param chrom: chromosome name (str)
        :param start: bin starting pos (int)
        :param end: bin end pos (int)
        :returns (str) sequence
        """
        return self.fasta.fetch(chrom, start, end)

    def __get_chr_lengths(self):
        """ returns dict with chromosome names and lengths
        :returns dictionary with chromosome name (str) and lengths(int)
        :rtype dictionary
        """
        lengths = self.fasta.lengths
        names = self.fasta.references
        
        return {name:length for name,length in zip(names, lengths)}
    
    def __get_chr_names(self):
        """extracts chromosome names from the bam file
        :returns list of chromosome names
        :rtype list 
        """
        return self.fasta.references
    
    
    def sliding_window(self, chrom, start, end):
        """simple sliding window implementation. init with 2*window_size bases.
        for each pos in the range (chrom start->end) return all chars from
        0->win_size and then remove one char from the beginning
        uses a deque for ease of insertion from right/deletion from left
        since pos increments with every loop and we pop once from left,
        start of the deque == pos
        
        :param chrom: (str) chromosome name
        :param start: (int) starting pos
        :param end: (int) end pos
        :returns generator object format: tuple(position, sequence)
        where sequence starts at position and is window_size long
        """
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
        """runs a sliding window from beginning to end of chromosome,
        prints sequences to output
        """
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
