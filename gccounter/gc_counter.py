'''
Created on Sep 22, 2017

@author: dgrewal
'''

import pysam
import argparse


class GCCounter( object):
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
    
    
    def get_gc_content(self, seq):
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
    
    parser.add_argument('reference')
    parser.add_argument('output')
    parser.add_argument('--chromosomes',
                        nargs='*')
    parser.add_argument('--window_size',
                        type=int,
                        default=1000)

    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    args = parse_args()

    GCCounter(args.reference, args.output, args.chromosomes, args.window_size).main()
