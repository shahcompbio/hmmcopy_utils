'''
Created on Aug 30, 2017

@author: dgrewal
'''
from subprocess import PIPE, Popen
import argparse

class bowtieIndex(object):
    def __init__(self, infile, outfile, reference, bwt_path):
        self.input = infile
        self.output = outfile
        self.bowtie = bwt_path
        self.reference = reference
     
    def cmd(self):
        cmd = [self.bowtie, self.reference,
               '-f', self.input, '-v',
               '0', '--quiet', '>', self.output]
     
        return cmd

    def main(self):
        cmd = self.cmd()

        cmd = " ".join(cmd)

        cmd = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
          
        stdout, stderr =  cmd.communicate()
        
        retc = cmd.returncode
        if retc != 0:
            raise Exception(stderr)
        
        return stdout

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('reference')
    parser.add_argument('--bowtie',
                        default="bowtie")

    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    args = parse_args()

    bowtieIndex(args.input, args.output, args.reference, args.bowtie).main()
