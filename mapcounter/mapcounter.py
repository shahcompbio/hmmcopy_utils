'''
Created on Aug 31, 2017

@author: dgrewal
'''
import pyBigWig as pybw
import argparse



class MapCounter(object):
    def __init__(self, infile, outfile, winsize):
        self.input = infile
        self.output = outfile
        self.window_size = winsize

    def gen_range(self, start, stop, step):
        current = start
        while current < stop:
            next_current = current + step
            if next_current < stop:
                yield (current, next_current)
            else:
                yield (current, stop)
            current = next_current

    def main(self):
        infile = pybw.open(self.input)
        output = open(self.output, 'w')
        
        chromosomes = infile.chroms()
        
        for chrom, len_chr in chromosomes.iteritems():
            #if window size > chromosome: raise err
            if self.window_size > len_chr:
                raise Exception("window size larger than chromosome")
            
            output.write("".join(map(str, ["fixedStep chrom=", chrom," start=1 step=",
                                    self.window_size, " span=", self.window_size])) + '\n')
            
            for start, stop in self.gen_range(1, len_chr, self.window_size):
                val = infile.values(chrom, start, stop)
                output.write("{0:.5f}".format(sum(val)/len(val)) + "\n")
#                 output.write(','.join(map(str,[start, sum(val), len(val)])) + '\n')

        infile.close()
        output.close()
        


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--window_size',
                        type=int,
                        default=1000)
    parser.add_argument('--seg',
                        default=False,
                        action="store_true")

    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    args = parse_args()

    MapCounter(args.input, args.output, args.window_size).main()