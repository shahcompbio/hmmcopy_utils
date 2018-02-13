'''
Created on Sep 29, 2017

@author: dgrewal
'''
from subprocess import Popen, PIPE
from correct_read_count import CorrectReadCount

from PyPDF2 import PdfFileMerger, PdfFileReader, PdfFileWriter

import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr

import os
import argparse
import math
import warnings

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def normalize_reads(df):
    if 'valid' not in df.columns.values:
        df['norm'] = float('nan')
    else:
        df['norm'] = df['reads']/np.median(df['reads'])
        df['norm'] = df['norm'].where(df['valid'] == True)
    return(df)


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


        self.mode = self.smoothing_function
        if self.mode == 'polynomial':
            self.mode += "_"+str(self.degree)


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
                c = CorrectReadCount(self.gc, self.map, inf, output, plot_out, sampid,
                                     mappability=args.mappability,
                                    smoothing_function=smoothing,
                                    polynomial_degree=degree)
                c.main()

        return plots, csvs

    def plot_pycorr(self, infile, outfile):
        pass

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

            skip=False
            if np.isinf(df["cor.map_hmmcopy"]).all() or np.isnan(df["cor.map_hmmcopy"]).all():
                print("%s sample skipped: hmmcopy map" %sampid)
                skip=True
            if np.isinf(df["cor.gc_hmmcopy"]).all() or np.isnan(df["cor.gc_hmmcopy"]).all():
                print("%s sample skipped: hmmcopy gc" %sampid)
                skip=True
            if np.isinf(df["cor.map_pycorr"]).all() or np.isnan(df["cor.map_pycorr"]).all():
                print("%s sample skipped: py map" %sampid)
                skip=True
            if np.isinf(df["cor.gc_pycorr"]).all() or np.isnan(df["cor.gc_pycorr"]).all():
                print("%s sample skipped: py gc" %sampid)
                skip=True
            if skip:
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
            hmm_data = normalize_reads(hmm_data)
            pycorr_data = pd.read_csv(pycorr[sampid])
            pycorr_data = normalize_reads(pycorr_data)

            pycorr_data = pycorr_data[pycorr_data['ideal'] == True]
            hmm_data = hmm_data[hmm_data["ideal"]==True]

            pycorr_data = pycorr_data.dropna()
            hmm_data = hmm_data.dropna()

            if pycorr_data.empty:
                continue

            if hmm_data.empty:
                continue

            #plot hmm
            # Four axes, returned as a 2-d array
            fig, axarr = plt.subplots(3, 2, sharex='col', sharey='row')
            fig.set_size_inches(10,15)


            axarr[0,0].scatter(pycorr_data["gc"], pycorr_data["reads"], alpha=0.1, s=4)
            axarr[0, 0].set_title('Uncorrected', fontsize=8)

            axarr[0,1].scatter(pycorr_data["map"], pycorr_data["reads"], alpha=0.1, s=4)
            axarr[0, 1].set_title('Uncorrected', fontsize=8)

            not_null = pycorr_data['cor.gc'].notnull()
            axarr[1,0].scatter(pycorr_data["gc"][not_null], pycorr_data["cor.gc"][not_null], alpha=0.1, s=4)
            axarr[1, 0].set_title('corrected ' + self.mode, fontsize=8)

            axarr[1,1].scatter(pycorr_data["map"][not_null], pycorr_data["cor.map"][not_null], alpha=0.1, s=4)
            axarr[1, 1].set_title('corrected ' + self.mode, fontsize=8)

            not_null = hmm_data['cor.gc'].notnull()
            axarr[2,0].scatter(hmm_data["gc"][not_null], hmm_data["cor.gc"][not_null], alpha=0.1, s=4)
            axarr[2, 0].set_title('corrected hmm', fontsize=8)

            axarr[2,1].scatter(hmm_data["map"][not_null], hmm_data["cor.map"][not_null], alpha=0.1, s=4)
            axarr[2, 1].set_title('corrected hmm', fontsize=8)

            axarr[2,0].set_xlabel("gc", fontsize=8)
            axarr[2,1].set_xlabel("map", fontsize=8)
            
            axarr[0,0].set_ylabel("reads", fontsize=8)
            axarr[1,0].set_ylabel(" normalized reads", fontsize=8)
            axarr[2,0].set_ylabel(" normalized reads", fontsize=8)

            plt.suptitle(sampid, fontsize=8)
            
            plots.savefig()
            plt.close()
        plots.close()

    def merge_by_cell(self, hmm_plots, py_plots, tempdir, output):
        samples = [v.replace('.wig','') for v in hmm_plots.keys()]

        for sample in samples:
            hmmplot = os.path.join(tempdir, sample+".wig_hmmcopy.pdf")
            pyplot = os.path.join(tempdir,sample+".wig_py.pdf")

            outpdf = os.path.join(tempdir,sample+"_merged.pdf")

            cmd = ['pdfjam', hmmplot, pyplot, '--nup 2x1', '--landscape', '--outfile', outpdf]

            self.run_cmd(cmd)

        merger = PdfFileMerger()

        for infile in samples:
            infile = os.path.join(tempdir,infile+"_merged.pdf")
            merger.append(open(infile, 'rb'))

        with open(output, 'wb') as fout:
            merger.write(fout)

    def main(self):
        hmm_plots, hmm_csv = self.run_hmmcopy()
        py_plots, py_csv = self.run_py_read_corr(self.smoothing_function, self.degree)

        self.merge_by_cell(hmm_plots, py_plots, self.tempdir, os.path.join(self.outdir, 'corr.pdf'))

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
    
