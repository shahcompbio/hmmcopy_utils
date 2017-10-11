'''
Created on Sep 22, 2017

@author: dgrewal
'''

import pysam
import argparse


class GCCounter( object):
    """
    calculates gc scores for segments of specified size
    """
    def __init__(self, ref_file, output, chromosomes, window_size):
        """
        :param ref_file: (str) reference fasta file
        :param output: (str) output wig file path
        :param chromosomes: list of target chromosomes (list of strings)
        :param window_size: (int) window sizre for binning genome
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
    
    
    def get_gc_content(self, seq):
        """calculates gc content from a list
        :param seq: (str) sequence
        :returns (float) gc_content for the sequence
        """
        length = 0
        gc = 0
        n =0 
         
        for val in seq:
            if val == 'G' or val == 'C' or val == 'g' or val == 'c':
                gc +=1
                length +=1  
            elif val == 'A' or  val == 'T' or val == 'a' or  val  == 't':
                length +=1
            elif val == '>':
                break
            else:
                n +=1
                length +=1
                break


        invalid = 1 if n > 0 else 0
        gc_content = -1 if invalid else float(gc)/(length-n)
#         n_content = float(n)/length

        gc_content= round(gc_content, 6)

        return gc_content
    
    def main(self):
        """calculates gc content for bins
        """
        with open(self.output, 'w') as outfile:
            for chromosome in self.chromosomes:
                outfile.write('fixedStep chrom=%s start=1 step=%s span=%s \n' %(chromosome, self.window_size, self.window_size))
                
                if chromosome not in self.chr_lengths:
                    raise Exception("Invalid chromosome")
            
                pos = 0
                while pos < self.chr_lengths[chromosome]:
                    seq = self.__get_sequence(chromosome, pos, pos + self.window_size)

                    gc = self.get_gc_content(seq)
    
                    outfile.write(str(gc) + '\n')
                    
                    pos += self.window_size



def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('reference',
                        help="path to the reference fasta file. The reference"\
                       " must be indexed with bowtie-build and samtools index")
                        
    parser.add_argument('output',
                        help="path to the output file")

    parser.add_argument('--chromosomes',
                        nargs='*',
                        help="specify target chromosomes")

    parser.add_argument('--window_size',
                        type=int,
                        default=1000,
                        help="specify window size.")

    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    args = parse_args()

    GCCounter(args.reference, args.output, args.chromosomes, args.window_size).main()
