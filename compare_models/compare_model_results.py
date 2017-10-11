'''
Created on Sep 29, 2017

@author: dgrewal
'''
from subprocess import Popen, PIPE
from correct_read_count import CorrectReadCount
from PyPDF2 import PdfFileMerger

import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr

import os
import argparse
import math


import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class Compare(object):
    def __init__(self, infiles, gc, mapp, tempdir, outdir, mappability=0.9, smoothing_function='lowess', degree=2):
        self.infiles = infiles

        self.gc = gc
        self.map = mapp

        self.tempdir = tempdir
        self.outdir = outdir
        
        self.mappability = mappability

        if not os.path.exists(tempdir):
            os.makedirs(tempdir)

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        self.smoothing_function = smoothing_function
        self.degree = degree

    @staticmethod
    def merge_pdf(infiles, outfile):
        merger = PdfFileMerger()
    
        for _,infile in infiles.iteritems():
            merger.append(open(infile, 'rb'))
    
        with open(outfile, 'wb') as fout:
            merger.write(fout)

    def run_cmd(self, cmd):
        """
        run a command with subprocess and return outputs
        """
        
        cmd = ' '.join(cmd)
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)

        _, stderr = proc.communicate()
        retc = proc.returncode

        if 'fatal' in stderr:
            raise Exception("The command %s failed with the "
                            "following error: %s" % (cmd, stderr))
        if retc != 0:
            raise Exception("The command %s failed with the "
                            "following error: %s" % (cmd, stderr))

    def run_hmmcopy(self):
        """
        runs hmm copy equiv script included
        """
        plots = {}
        csvs = {}
        
        for inf in self.infiles:
            sampid = os.path.basename(inf)


            output = os.path.join(self.tempdir, sampid+"_hmmcopy.csv")
            plot_out = os.path.join(self.tempdir, sampid+"_hmmcopy.pdf")
            
            plots[sampid] =  plot_out
            csvs[sampid] = output

            if not os.path.exists(plot_out) or not os.path.exists(output): 
                cmd = ["Rscript", "hmmcopy_corr.R", inf, self.map, self.gc, output, plot_out,
                       str(self.mappability), sampid]
                
                self.run_cmd(cmd)

        return plots, csvs

    def run_py_read_corr(self, smoothing, degree):
        plots = {}
        csvs = {}
        
        for inf in self.infiles:

            sampid = os.path.basename(inf)
            output = os.path.join(self.tempdir, sampid+"_py.csv")
            plot_out = os.path.join(self.tempdir, sampid+"_py.pdf")

            plots[sampid] =  plot_out
            csvs[sampid] = output

            if not os.path.exists(plot_out) or not os.path.exists(output):
                c = CorrectReadCount(self.gc, self.map, inf, output, plot_out,
                                     sampid, mappability=args.mappability,
                                    smoothing_function=smoothing,
                                    polynomial_degree=degree)
                c.main()

        return plots, csvs

    def plot_df(self, pearson_data, mean_data):
        plots = os.path.join(self.outdir, 'plots.pdf')
        plots = PdfPages(plots)
        
        #plot correlation
        plt.Figure()
        gc = list(zip(*pearson_data))[0]
        mapp = list(zip(*pearson_data))[1]
        plt.plot(gc, label='gc')
        plt.plot(mapp, label='mapp')
        plt.title("pearson")
        plt.legend()
        plots.savefig()
        plt.close()

        #plot min
        plt.Figure()
        gc = list(zip(*mean_data))[0]
        mapp = list(zip(*mean_data))[1]
        plt.plot(gc, label='gc')
        plt.plot(mapp, label='mapp')
        plt.title("mean diff")
        plt.legend()
        plt.ylim(0,1)
        plots.savefig()

        plt.close()
        plots.close()

    def plot_metrics(self, hmm, pycorr):

        pearson_data = []
        mean_diff = []
        
        for sampid, hmm_data in hmm.iteritems():
            hmm_data = pd.read_csv(hmm_data)
            pycorr_data = pd.read_csv(pycorr[sampid])
            
            df = pd.merge(hmm_data, pycorr_data,
                          on=['chromosome','start','end', 'reads', 'gc', 'map'],
                          suffixes=["_hmmcopy", "_pycorr"])
            
            df = df.dropna() 
            
            if df.empty:
                continue           
            
            gc_corr =  pearsonr(df["cor.gc_hmmcopy"], df["cor.gc_pycorr"])[0]
            map_corr =  pearsonr(df["cor.map_hmmcopy"], df["cor.map_pycorr"])[0]

            pearson_data.append((gc_corr, map_corr))

            if np.isinf(df["cor.map_hmmcopy"]).all() or np.isnan(df["cor.map_hmmcopy"]).all():
                continue
            if np.isinf(df["cor.gc_hmmcopy"]).all() or np.isnan(df["cor.gc_hmmcopy"]).all():
                continue

            gc = [abs(a-b) for a,b in zip(df["cor.gc_hmmcopy"], df["cor.gc_pycorr"])]
            mapp = [abs(np.subtract(a,b)) for a,b in zip(df["cor.map_hmmcopy"], df["cor.map_pycorr"])]
            mean_diff.append((np.mean(gc), np.mean(mapp)))

        self.plot_df(pearson_data, mean_diff)


    def boxplots(self, hmm, pycorr):
        plots = os.path.join(self.outdir, 'boxplots.pdf')
        plots = PdfPages(plots)
        

        for sampid, hmm_data in hmm.iteritems():
            hmm_data = pd.read_csv(hmm_data)
            pycorr_data = pd.read_csv(pycorr[sampid])
            
            df = pd.merge(hmm_data, pycorr_data,
                          on=['chromosome','start','end', 'reads', 'gc', 'map'],
                          suffixes=["_hmmcopy", "_pycorr"])
            
            df = df.dropna() 
            
            if df.empty:
                continue           

            gc = [abs(np.subtract(a,b)) for a,b in zip(df["cor.gc_hmmcopy"], df["cor.gc_pycorr"]) if not math.isnan(np.subtract(a,b))]
            mapp = [abs(np.subtract(a,b)) for a,b in zip(df["cor.map_hmmcopy"], df["cor.map_pycorr"]) if not math.isnan(np.subtract(a,b))]
           
            if not gc:
                continue
            if not mapp:
                continue

            plt.Figure()
            plt.boxplot([gc,mapp], showfliers=False)
            
            plt.xticks( [1,2], ['gc','mapp'], rotation='vertical')
            
            plt.title(sampid)

            plots.savefig()
            plt.close()
        plots.close()

    def plot_scatter(self, hmm, pycorr):

        plots = os.path.join(self.outdir, 'scatterplots.pdf')
        plots = PdfPages(plots)
        

        for sampid, hmm_data in hmm.iteritems():
            hmm_data = pd.read_csv(hmm_data)
            pycorr_data = pd.read_csv(pycorr[sampid])
            
            df = pd.merge(hmm_data, pycorr_data,
                          on=['chromosome','start','end', 'reads', 'gc', 'map'],
                          suffixes=["_hmmcopy", "_pycorr"])
            
            df = df.dropna() 
            
            if df.empty:
                continue           


            plt.Figure()
            #plot hmm
            # Four axes, returned as a 2-d array
            _, axarr = plt.subplots(2, 2, sharex=True, sharey=True)
            
            
            axarr[0,0].scatter(df["reads"], df["cor.gc_pycorr"], alpha=0.01, s=4)
            axarr[0, 0].set_title('reads vs gc corrected py')
            axarr[0,1].scatter(df["reads"], df["cor.gc_hmmcopy"], alpha=0.01, s=4)
            axarr[0, 1].set_title('reads vs gc corrected hmmcopy')
            
            axarr[1,0].scatter(df["reads"], df["cor.map_pycorr"], alpha=0.01, s=4)
            axarr[1, 0].set_title('reads vs map corrected py')
            axarr[1,1].scatter(df["reads"], df["cor.map_hmmcopy"], alpha=0.01, s=4)
            axarr[1, 1].set_title('reads vs map corrected hmmcopy')
            
            axarr[0,0].set_ylim((0,2.5))
            axarr[0,0].set_ylabel("gc")
            axarr[1,0].set_ylabel("map")
            
            axarr[1,0].set_xlabel("reads")
            axarr[1,1].set_xlabel("reads")
            
            plt.suptitle(sampid)
            
            plots.savefig()
            plt.close()
        plots.close()

    def main(self):
        hmm_plots, hmm_csv = self.run_hmmcopy()
        py_plots, py_csv = self.run_py_read_corr(self.smoothing_function, self.degree)
        
        self.merge_pdf(hmm_plots, os.path.join(self.outdir, 'hmm.pdf'))
        self.merge_pdf(py_plots, os.path.join(self.outdir, 'py.pdf'))

        #some comparison metrics
        self.plot_metrics(hmm_csv, py_csv)

        #boxplots of diff.
        self.boxplots(hmm_csv, py_csv)
        self.plot_scatter(hmm_csv, py_csv)


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--infiles',
                        nargs = "*")

    parser.add_argument('--gc')
    
    parser.add_argument('--map')

    parser.add_argument('--tempdir')
    parser.add_argument('--outdir')
    parser.add_argument('--mappability',
                        default=0.9)
    
    parser.add_argument('--smoothing_function',
                        default='lowess')
    
    parser.add_argument('--degree',
                        default=2,
                        type=int)


    
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':
    args = parse_args()
    
    comp = Compare(args.infiles, args.gc, args.map, args.tempdir, args.outdir, mappability=args.mappability,
                   smoothing_function=args.smoothing_function, degree=args.degree)
    comp.main()
    
