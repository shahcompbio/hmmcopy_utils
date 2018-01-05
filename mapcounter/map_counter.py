'''
Created on Aug 29, 2017

@author: dgrewal
'''
import argparse
import pysam
from collections import defaultdict


class CountMappability(object):
    """generates mappability bigwig file from the aligned bowtie output
    """

    def __init__(
            self, infile, outfile, reference, chromosomes, maxhits, winsize):
        """
        :param infile:input bowtie file
        :param infile:output mappability bigwig file
        :param maxhits: (int) max hits
        """
        self.window_size = winsize
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

        return {name: length for name, length in zip(names, lengths)}

    def get_bins(self, chrom):
        """divide the input region (start,stop) into step sized regions
        :param start: (int) start position
        :param stop: (int) stop position
        :param step: (int) window size
        """
        length = self.chr_lengths[chrom]
        step = self.window_size

        bins = []

        current = 1
        while current < length:
            next_current = current + step
            if next_current < length:
                bins.append((current, next_current))
            else:
                bins.append((current, length))
            current = next_current

        return bins

    def get_bin_val(self, chrom, pos):

        step = self.window_size
        start = (pos / step) * step
        stop = start + step

        start += 1
        stop += 1

        stop = min(self.chr_lengths[chrom], stop)

        return start, stop

    def parse_bowtie(self):
        """parse bowtie output, assumes the input file is sorted.
        """
        data = defaultdict(lambda: defaultdict(float))
        expected_pos = 0

        with open(self.input) as infile:
            for line in infile:
                line = line.strip()

                if line == "":
                    continue

                pos, _, _, _, _, _, hits = line.split()
                chrom, pos = pos.split(":")

                if chrom not in self.chromosomes:
                    continue

                pos = int(pos)

                bin_val = self.get_bin_val(chrom, pos)

                hits = float(hits)
                hits += 1
                value = 0.0
                if hits <= self.maxhits:
                    value = 1.0 / hits

                data[chrom][bin_val] += value
                expected_pos += 1

        return data

    def write_wig(self, data):

        with open(self.output, 'w') as outfile:

            for chrom, bin_data in data.iteritems():
                outfile.write("".join(map(str, ["fixedStep chrom=", chrom, " start=1 step=",
                                                self.window_size, " span=", self.window_size])) + '\n')

                for binval in self.get_bins(chrom):
                    binsize = binval[1] - binval[0]
                    value = float(bin_data[binval]) / binsize
                    outfile.write("{0:.5f}".format(value) + '\n')

    def main(self):
        data = self.parse_bowtie()

        self.write_wig(data)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('reference')
    parser.add_argument('--chromosomes',
                        nargs="*")
    parser.add_argument("--window_size",
                        type=int)
    parser.add_argument('--max_hits',
                        type=int,
                        default=4)

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()

    CountMappability(
        args.input,
        args.output,
        args.reference,
        args.chromosomes,
        args.max_hits,
        args.window_size).main()
