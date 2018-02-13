'''
Created on Sep 27, 2017


@author: dgrewal
'''
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.interpolate import interp1d
from scipy.stats.mstats import mquantiles
import numpy.polynomial.polynomial as poly

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class CorrectReadCount(object):
    """
    fit lowess/polynomial curve smoothing to reads-gc
    use fitted model to predict corrected gc and mappability
    values
    """

    def __init__(self, gc, mapp, wig, output, plot_out, sample_id, mappability=0.9,
                 smoothing_function='lowess',
                 polynomial_degree=2):
        self.mappability = mappability
        self.smoothing_function = smoothing_function
        self.polynomial_degree = polynomial_degree

        self.gc = gc
        self.mapp = mapp
        self.wig = wig
        self.output = output
        self.sampid = sample_id

        self.pdfout = PdfPages(plot_out)

        self.mode = self.smoothing_function
        if self.mode == 'polynomial':
            self.mode += "_"+str(self.polynomial_degree)


    def plot(self, scatter_data, line_data, title, alpha=0.05):

        title = self.sampid + ' '+ title

        plt.Figure(figsize=(4,6))
        plt.scatter(scatter_data[1], scatter_data[0], alpha=alpha)
        plt.plot(line_data[0], line_data[1])
        plt.title(title)
        plt.ylabel("gc")
        plt.xlabel("reads")
        
        self.pdfout.savefig()
        plt.close()

    def read_wig(self, infile, counts=False):
        """read wiggle files

        :param infile: input wiggle file
        :param counts: set to true if infile wiggle has integer values
        """

        data = []

        with open(infile) as wig:
            for line in wig:
                line = line.strip()

                if line.startswith('fixedStep'):
                    line = line.strip().split()

                    chrom = line[1].split('=')[1]
                    winsize = int(line[3].split('=')[1])
                    start = int(line[2].split('=')[1])

                    bin_start = 0 if start < winsize else start / winsize
                else:
                    value = int(line) if counts else float(line)
                    data.append((chrom, (bin_start * winsize) + 1,
                                 (bin_start + 1) * winsize, winsize, value))
                    bin_start += 1

        return data

    def create_dataframe(self, reads, mapp, gc):
        """merge data from reads, mappability and gc wig files
        into pandas dataframe

        :param reads: list of tuples, formatted as [(chromosome,
                    start, end, count), ]
        :param mapp: list of tuples, formatted as [(chromosome,
                    start, end, mappability value), ]
        :param reads: list of tuples, formatted as [(chromosome,
                    start, end, gc content), ]
        """
        err_str = 'please ensure that reads, mappability and '\
            'gc wig files have the same sort order'

        data = []
        for read_v, mapp_v, gc_v in zip(reads, mapp, gc):
            assert read_v[0] == mapp_v[0] == gc_v[0], err_str
            assert read_v[1] == mapp_v[1] == gc_v[1], err_str
            assert read_v[2] == mapp_v[2] == gc_v[2], err_str
            assert read_v[3] == mapp_v[3] == gc_v[3], err_str

            data.append((read_v[0], read_v[1], read_v[2], read_v[3], gc_v[4],
                         mapp_v[4], read_v[4],))

        labels = ['chromosome', 'start', 'end', 'width', 'gc', 'map', 'reads']
        data = pd.DataFrame(data, columns=labels)

        return data

    def valid(self, df):
        """adds valid column (calls with atleast one reads and non negative gc)

        :params df: pandas dataframe
        """

        df.loc[:, "valid"] = True

        df.loc[(df["reads"] <= 0) | (df['gc'] < 0), "valid"] = False

        return df

    def ideal(self, df):
        """adds ideal column

        :params df: pandas dataframe
        """
        df.loc[:, "ideal"] = True

        valid_reads = df[df["valid"]]["reads"]
        valid_gc = df[df["valid"]]["gc"]

        routlier = 0.01
        doutlier = 0.001

        range_l, range_h = mquantiles(valid_reads, prob=[0, 1 - routlier],
                                      alphap=1, betap=1)
        domain_l, domain_h = mquantiles(valid_gc, prob=[doutlier, 1 - doutlier],
                                        alphap=1, betap=1)

        df.loc[(df["valid"] == False) |
               (df["map"] < self.mappability) |
               (df["reads"] <= range_l) |
               (df["reads"] > range_h) |
               (df["gc"] < domain_l) |
               (df["gc"] > domain_h),
               "ideal"] = False

        return df

    def get_lowess_fit_gc(self, x, y, x_):
        """fits lowess curve to [x,y]

        :param x: numpy array
        :param y: numpy array
        """

        #sm lowess takes endog,exog as input, which translates to y,x
        #returns sorted x values and corresponding smoothed y values
        lowess = sm.nonparametric.lowess(y, x, frac=.03)
        # unpack the lowess smoothed points to their values
        lowess_x = list(zip(*lowess))[0]
        lowess_y = list(zip(*lowess))[1]

        #plot start
        scatter_data = (x,y)
        line_data = (lowess_y, lowess_x)
        self.plot(scatter_data, line_data, self.mode)
        #plot end

        f = interp1d(lowess_x, lowess_y, bounds_error=False)

        i = np.arange(0, 1, 0.001)

        preds = f(i)

        final_lowess = sm.nonparametric.lowess(preds, i, frac=0.3)
        lowess_x = list(zip(*final_lowess))[0]
        lowess_y = list(zip(*final_lowess))[1]


#         #plot start
#         scatter_data = (preds, i)
#         line_data = (lowess_y, lowess_x)
#         self.plot(scatter_data, line_data, 'pycorrect second pass')
#         #plot end

        f = interp1d(lowess_x, lowess_y, bounds_error=False)


        return f(x_)

    def get_poly_fit(self, x, y, x_):
        """fits polynomial curve to [x,y] with degree=self.polynomial_degree

        :param x: numpy array
        :param y: numpy array
        """

        coefs = poly.polyfit(x,y, self.polynomial_degree)
        c =  poly.polyval(x_, coefs)


        c_p = poly.polyval(x, coefs)
        #plot start
        scatter_data = (x,y)

        mapping ={y:x for x,y in zip(c_p, x)}
        c_p = [mapping[v] for v in sorted(x)]
        line_data = (c_p, sorted(x))
        self.plot(scatter_data, line_data, self.mode)
         #plot end

        return c

    def get_lowess_fit_mapp(self, x, y, x_):
        """fits lowess curve to [x,y]

        :param x: numpy array
        :param y: numpy array
        """
        delta=0.01*(max(x) - min(x))
        #x:map y:gc
        lowess = sm.nonparametric.lowess(y,x, frac=.66, delta = delta)
        # unpack the lowess smoothed points to their values
        lowess_x = list(zip(*lowess))[0]
        lowess_y = list(zip(*lowess))[1]

        f = interp1d(lowess_x, lowess_y, bounds_error=False)

        return f(x_)

    def correct_gc(self, df):
        """calculate corrected gc values with the specified model

        :param df: pandas dataframe
        """

        # only keep ideal val
        ideal_df = df.loc[df["ideal"] == True]

        gc = np.array(ideal_df["gc"])
        reads = np.array(ideal_df["reads"])

        if self.smoothing_function == 'polynomial':
            preds = self.get_poly_fit(gc, reads, df["gc"])
        else:
            preds = self.get_lowess_fit_gc(gc, reads, df["gc"])

        df.loc[:, 'cor.gc'] = df["reads"] / preds

        return df

    def correct_mapp(self, df):
        """calculate corrected map values with the specified model

        :param df: pandas dataframe
        """

        valid_gc = df[df["valid"]]["cor.gc"].dropna()

        coutlier = 0.01
        range_h = mquantiles(valid_gc, prob=1 - coutlier,
                             alphap=1, betap=1)
        range_h = range_h[0]

        ideal_df = df.loc[(df["cor.gc"] < range_h)]

        preds = self.get_lowess_fit_mapp(ideal_df["map"], ideal_df["cor.gc"], df["map"])

        df.loc[:, "cor.map"] = df["cor.gc"] / preds
        df.loc[:, "copy"] = np.log2(df['cor.map'])

        return df

    def write(self, df):
        """write results to the output file

        :param df: pandas dataframe
        """

        df.to_csv(self.output, index=False, sep=',', na_rep="NA")

    def main(self):
        gc = self.read_wig(self.gc)
        mapp = self.read_wig(self.mapp)
        reads = self.read_wig(self.wig, counts=True)

        df = self.create_dataframe(reads, mapp, gc)

        df = self.valid(df)
        df = self.ideal(df)

        df = self.correct_gc(df)
        df = self.correct_mapp(df)

        self.write(df)

        self.pdfout.close()


def parse_args():
    """
    parses command line arguments
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('gc',
                        help='path to the gc wig file'
                        )

    parser.add_argument('map',
                        help='path to the mappability wig file'
                        )

    parser.add_argument('reads',
                        help='path to the read-counts wig file'
                        )

    parser.add_argument('output',
                        help='path to the output csv file'
                        )

    parser.add_argument('--smoothing_function',
                        default='lowess',
                        choices=['lowess', 'polynomial'],
                        help='specify the curve fitting algorithm')

    parser.add_argument('--mappability',
                        default=0.9,
                        type=float,
                        help='specify mappability threshold')

    parser.add_argument('--degree',
                        default=2,
                        type=int,
                        help='specify degree for the polynomial fit, only '\
                        'used if the smoothing function is set to polynomial')


    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()

    corr = CorrectReadCount(args.gc, args.map, args.reads, args.output,
                            mappability=args.mappability,
                            smoothing_function=args.smoothing_function,
                            polynomial_degree=args.degree
                            )

    corr.main()
