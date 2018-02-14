'''
Created on Oct 10, 2017

@author: dgrewal
'''
import os
import pysam
import argparse


class ReadCounter(object):
    """
    calculate reads per bin from the input bam file
    """

    def __init__(self, bam, output, window_size, chromosomes, mapq, seg):
        self.bam = bam

        self.output = output

        self.window_size = window_size

        if chromosomes:
            self.chromosomes = chromosomes
        else:
            self.chromosomes = self.__get_chr_names()

        self.bam = self.__get_bam_reader()
        self.chr_lengths = self.__get_chr_lengths()

        self.mapq_threshold = mapq

        self.seg = seg

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        #clean up output if there are any exceptions
        if exc_type and os.path.exists(self.output):
            os.remove(self.output)


    def __get_chr_lengths(self):
        """ returns dict with chromosome names and lengths
        :returns dictionary with chromosome name (str) and lengths(int)
        :rtype dictionary
        """

        names = self.bam.references
        lengths = self.bam.lengths
        return {name: length for name, length in zip(names, lengths)}

    def __get_bam_reader(self):
        """returns pysam bam object
        :returns pysam bam object
        """
        return pysam.AlignmentFile(self.bam, 'rb')

    def __get_chr_names(self):
        """extracts chromosome names from the bam file
        :returns list of chromosome names
        :rtype list 
        """
        return self.bam.references

    def __fetch(self, chrom, start, end):
        """returns iterator over reads in the specified region
        :param chrom: chromosome name (str)
        :param start: bin starting pos (int)
        :param end: bin end pos (int)
        :returns iterator over reads
        """
        return self.bam.fetch(chrom, start, end)

    def filter(self, pileupobj):
        """remove low mapping quality reads and duplicates
        :param pileupobj: pysam read object
        :returns boolean: false if the read passes filters.
        :rtype boolean
        """

        if not pileupobj.is_duplicate and \
                pileupobj.mapping_quality >= self.mapq_threshold:
            return False

        return True

    def write_header(self, chrom, outfile):
        """writes headers, single header if seg format,
        one header per chromosome otherwise.
        :param chrom: chromosome name
        :param outfile: output file object.
        """
        if self.seg:
            outfile.write("data\tchr\tstart\tend\tcount\n")
        else:
            outstr = "fixedStep chrom=%s start=1 step=%s span=%s\n"\
                % (chrom, self.window_size, self.window_size)
            outfile.write(outstr)

    def write(self, chrom, start, stop, count, outfile):
        """writes bin and counts to the output file.
        supports seg and wig formats
        :param chrom: chromosome name
        :param start: bin start
        :param stop: bin stop
        :param count: no of reads in the bin
        :param outfile: output file object.
        """
        if self.seg:
            outstr = '\t'.join(
                map(str, ['reads', chrom, start, stop, count])) + '\n'
            outfile.write(outstr)
        else:
            outfile.write(str(count) + '\n')

    def get_data(self, data, chrom, outfile):
        """iterates over reads, calculates counts and writes to output
        :param data: pysam iterator over reads
        :param chrom: str: chromosome name
        :param outfile: output file object
        """
        reflen = self.chr_lengths[chrom]

        count = 0
        start = 0
        end = start + self.window_size
        for pileupobj in data:
            while pileupobj.pos > end:
                self.write(chrom, start, end, count, outfile)
                count = 0
                start += self.window_size
                end = min(start + self.window_size, reflen)

            if not self.filter(pileupobj):
                count += 1

        while end < reflen:
            self.write(chrom, start, end, 0, outfile)
            start += self.window_size
            end = min(start + self.window_size, reflen)
        self.write(chrom, start, end, 0, outfile)



    def main(self):
        """for each chromosome, iterate over all reads. use starting position
        of the read to calculate read counts per bin (no double counting).
        """
        with open(self.output, 'w') as outfile:
            if self.seg:
                self.write_header(None, outfile)

            for chrom in self.chromosomes:
                if not self.seg:
                    self.write_header(chrom, outfile)

                reflen = self.chr_lengths[chrom]

                # get read iterator for the full chromosome
                # code assumes the iterator is sorted.
                data = self.__fetch(chrom, 0, reflen)

                self.get_data(data, chrom, outfile)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('bam',
                        help='specify the path to the input bam file')

    parser.add_argument('output',
                        help='specify path to the output file')
    parser.add_argument('--chromosomes',
                        nargs='*',
                        default=map(str, range(1, 23)) + ['X', 'Y'],
                        help='specify target chromosomes'
                        )
    parser.add_argument('-w', '--window_size',
                        type=int,
                        default=1000,
                        help='specify bin size')
    parser.add_argument('-m', '--mapping_quality_threshold',
                        type=int,
                        default=0,
                        help='threshold for the mapping quality, reads '\
                        'with quality lower than threshold will be ignored')

    parser.add_argument('--seg',
                        default=False,
                        action='store_true',
                        help='write the output in seg format')

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = parse_args()
    with ReadCounter(args.bam, args.output, args.window_size,
                     args.chromosomes, args.mapping_quality_threshold,
                     args.seg) as rcount:
        rcount.main()
