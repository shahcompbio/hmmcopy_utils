'''
Created on Aug 29, 2017

@author: dgrewal
'''
import argparse
import pyBigWig as pybw
import pysam
import numpy as np
import gzip

class readToMap(object):
    """generates mappability bigwig file from the aligned bowtie output
    """

    def __init__(self, infile, outfile, reference, chromosomes, maxhits):
        """
        :param infile:input bowtie file
        :param infile:output mappability bigwig file
        :param maxhits: (int) max hits
        """
        self.input = infile
        self.output = outfile
        self.maxhits = maxhits
        self.chr_lengths = self.__get_chr_lengths(reference)
        self.chromosomes = chromosomes

    def __get_chr_lengths(self, reference):
        """ returns dict with chromosome names and lengths
        :returns dictionary with chromosome name (str) and lengths(int)
        :rtype dictionary
        """
        fasta = pysam.FastaFile(reference)

        lengths = fasta.lengths
        names = fasta.references

        return {name:length for name,length in zip(names, lengths)}

    def add_header_bigwig(self, bigwig):


        header = [(chrom, self.chr_lengths[chrom]) for chrom in self.chromosomes]

        bigwig.addHeader(header)

    def parse_bowtie(self):
        """parse bowtie output, assumes the input file is sorted.
        """

        last_chrom = None
        expected_pos = 0

        values = []
        with gzip.open(self.input) as infile:
            for line in infile:
                line = line.strip()

                if line == "":
                    continue

                pos, _, _, _, _, _, hits = line.split()
                chrom, pos = pos.split(":")

                hits = float(hits)
                pos = int(pos)
                
                if not last_chrom:
                    last_chrom = chrom

                if not chrom == last_chrom:
                    yield last_chrom, values
                    last_chrom = chrom
                    expected_pos = 0
                    values = []

                if expected_pos <= pos:
                    if values:
                        yield chrom,values
                        values = []

                    diff = pos - expected_pos
                    for _ in range(diff/1000000):
                        yield chrom, np.zeros(1000000)
                    yield chrom, np.zeros(diff%1000000)
                    expected_pos += diff

                hits += 1
                value = 0.0
                if hits <= self.maxhits:
                    value = 1.0 / hits

                values.append(value)
                expected_pos += 1

                if len(values) >= 1000000:
                    yield chrom, values
                    values = []

            yield chrom, values

    def main(self):
        output = pybw.open(self.output, "w")

        self.add_header_bigwig(output)

        vals = self.parse_bowtie()

        last_chrom = None
        pos = 1

        for chrom, values in vals:
            if not last_chrom:
                last_chrom = chrom

            #add header if chrom changes
            if last_chrom != chrom:
                pos = 1
                last_chrom = chrom

            output.addEntries(chrom, pos, values=values, span=1, step=1)
            pos += len(values)
        output.close()

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('reference')
    parser.add_argument('--chromosomes',
                        nargs="*")
    parser.add_argument('--max_hits',
                        type=int,
                        default=4)

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()

    readToMap(args.input, args.output, args.reference, args.chromosomes, args.max_hits).main()
