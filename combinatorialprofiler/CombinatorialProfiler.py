#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
import os.path
import subprocess
import csv
import json
import logging
import time
import pickle

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('PDF')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import gaussian_kde

import Bio.Seq
import Bio.Alphabet

from .readcounter import PyHammingReadCounter, PySeqlevReadCounter, PyExperiment, DSIS
from . import version

def formatTime(seconds):
    if seconds < 60:
        return "%i seconds" % round(seconds)
    elif seconds < 3600:
        return "%.1f minutes" % (seconds / 60)
    else:
        return "%.1f hours" % (seconds / 3600)

class TimeLogger:
    def __init__(self, startmsg, stopmsg=None, level=logging.INFO):
        self._startmsg = startmsg[:1].upper() + startmsg[1:]
        if stopmsg:
            self._stopmsg = stopmsg[:1].upper() + stopmsg[1:]
        else:
            self._stopmsg = "Finished %s" % startmsg
        self._level = level
        self._appendstopmsg = ""

    def __enter__(self):
        self._starttime = time.monotonic()
        logging.log(self._level, self._startmsg)

    def __exit__(self, *args):
        stop = time.monotonic()
        logging.log(self._level, "%s after %s%s" % (self._stopmsg, formatTime(stop - self._starttime), self._appendstopmsg))

    def setLevel(level):
        self._level = level

    def setStopMsg(msg):
        self._stopmsg = msg

    def appendToStopMsg(msg):
        self._appendstopmsg = msg

def normalizeCounts(df, sortedcells):
    categories = (df['experiment'].cat.categories.sort_values(), df['barcode_fw'].cat.categories.sort_values(), df['barcode_rev'].cat.categories.sort_values())
    df = df.set_index(['experiment','barcode_fw', 'barcode_rev','sequence'])
    df['normalized_counts'] = (df['counts'] / df.groupby(level=['barcode_fw', 'barcode_rev'])['counts'].transform(sum)).unstack(['experiment','sequence']).mul(sortedcells, axis=0).stack(['experiment','sequence']).reorder_levels(df.index.names).dropna()
    df = df.reset_index()
    df['experiment'] = df['experiment'].astype("category", categories=categories[0])
    df['barcode_fw'] = df['barcode_fw'].astype("category", categories=categories[1])
    df['barcode_rev'] = df['barcode_rev'].astype("category", categories=categories[2])
    return df

class DSISpec(object):
    pass

def getDSISpec(e):
    dspec = DSISpec()
    if isinstance(e, PyExperiment):
        dsi = e.dsi
        if len(e.named_inserts):
            dspec.seqcol = 'named_insert'
        else:
            dspec.seqcol = 'sequence'
    elif isinstance(e, DSIS):
        dsi = e
    else:
        raise TypeError("Unsupported type")

    if dsi == DSIS.noDSI:
        return None
    elif dsi == DSIS.forward:
        dspec.groupby = 'barcode_rev'
        dspec.dsicol = 'barcode_fw'
    elif dsi == DSIS.reverse:
        dspec.groupby = 'barcode_fw'
        dspec.dsicol = 'barcode_rev'

    return dspec

def getDSI(df, dspec):
    statsfuns = ['min', 'max', 'mean', 'median', 'std', 'sum']

    fractionvals = pd.Series(range(1, df[dspec.dsicol].cat.categories.size + 1), index=sorted(df[dspec.dsicol].cat.categories))

    df['normalized_counts_cells'] = df.set_index(dspec.dsicol, append=True)['normalized_counts'].mul(fractionvals, level=dspec.dsicol).reset_index(level=dspec.dsicol, drop=True)

    groupbyl = ['experiment', dspec.groupby]
    if 'named_insert' in df.columns:
        groupbyl.append('named_insert')

    g = df.groupby(groupbyl + ['translation','sequence'])
    byseq = (g['normalized_counts_cells'].sum() / g['normalized_counts'].sum()).dropna()
    byseq.name = 'dsi'
    byseq_stats = g.agg({'counts': statsfuns, 'normalized_counts': statsfuns})
    byseq_stats.columns = ['_'.join(col) for col in byseq_stats.columns.values]
    byseq_stats['nfractions'] = g.agg('size')

    byaa_median = byseq.groupby(level=groupbyl + ['translation']).median().dropna()
    byaa_median.name = 'median_dsi'

    g = df.groupby(groupbyl + ['translation'])
    byaa_pooled = (g['normalized_counts_cells'].sum() / g['normalized_counts'].sum()).dropna()
    byaa_pooled.name = 'pooled_dsi'

    byaa_stats = g.agg({'counts': statsfuns, 'normalized_counts': statsfuns})
    byaa_stats.columns = ['_'.join(col) for col in byaa_stats.columns.values]
    byaa_stats['nfractions'] = g[dspec.dsicol].nunique()

    return (pd.concat((byaa_median, byaa_pooled, byaa_stats), axis=1).reset_index().dropna(subset=('median_dsi', 'pooled_dsi')), pd.concat((byseq, byseq_stats), axis=1).reset_index().dropna(subset=['dsi']))

class PyExperimentJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if type(o) != PyExperiment:
            return super().default(o)
        else:
            return {o.name : o.toDict()}

def exec_with_logging(args, pname, out=None, err=None):
    if out:
        outf = open(out, 'w')
    else:
        outf = sys.stdout
    if err:
        errf = open(err, 'w')
    else:
        errf = sys.stderr

    with TimeLogger("Starting %s" % pname, "%s finished" % pname) as tl:
        logging.debug(" ".join(args))
        ret = subprocess.call(args, stdout=outf, stderr=errf)
        if ret:
            tl.appendToStopMsg(" with returncode %d" % ret)
            tl.setLevel(logging.ERROR)
        return ret

def plot_profiles(df, groupby, dspec, filename):
    labels = sorted(df[dspec.dsicol].cat.categories)
    integer_map = dict([(val, i) for i, val in enumerate(labels)])
    with PdfPages(filename) as pdf:
        for (code, seq), group in df.groupby(groupby):
            fig = plt.figure(figsize=(5,3))

            xvals = group[dspec.dsicol].map(integer_map)
            xorder = xvals.argsort()
            plot = plt.plot(xvals.iloc[xorder], group['normalized_counts'].iloc[xorder], 'k.-')
            plt.xlim(0, len(labels) - 1)
            plt.ylim(ymin=0)
            plt.xticks(range(len(labels)), labels)
            plt.title("%s %s" % (code, seq))
            plt.ylabel("normalized counts")
            pdf.savefig(bbox_inches='tight')
            plt.close()

def plot_histograms(df, dspec, filename, quantile=1):
    with PdfPages(filename) as pdf:
        for code, group in df.groupby(dspec.groupby):
            counts = group['counts_sum'][group['counts_sum'] <= group['counts_sum'].quantile(quantile)]

            fig = plt.figure(figsize=(5,3))
            limits = (0, counts.max())
            nbins = min(100, int(counts.size / 10))
            hist = plt.hist(counts, nbins, range=limits, normed=True, color="#000000", alpha=0.66, edgecolor='none')
            plt.xlim(limits)
            plt.xlabel("Total read count")
            plt.ylabel("Frequency")

            smoothedbins = np.linspace(limits[0], limits[1], nbins * 10)
            kde = gaussian_kde(counts)
            kde.set_bandwidth(kde.factor * 0.75)
            y = kde.evaluate(smoothedbins)
            plt.plot(smoothedbins, y, 'k-')

            plt.title(code)
            pdf.savefig(bbox_inches='tight')
            plt.close()

def plot_correlations(df, dspec, limits, filename):
    with PdfPages(filename) as pdf:
        for code, group in df.groupby(dspec.groupby):
            c = group['median_dsi'].corr(group['pooled_dsi'], method='spearman')

            fig = plt.figure(figsize=(7,7))

            ax = fig.add_subplot(1,1,1)
            ax.scatter(group['median_dsi'], group['pooled_dsi'], s=100, c="#000000", alpha=1/3, marker='.', edgecolor='none')

            ax.set_xlabel("median DSI")
            ax.set_ylabel("pooled DSI")
            ax.text(0.1, 0.9, "Spearman's $r = %.3g$" % c, transform=ax.transAxes)
            ax.set_aspect('equal', adjustable='box', anchor='C')

            nbins = (limits[1] - limits[0] + 1) * 5
            binwidth = (limits[1] - limits[0] + 1) / nbins
            divider = make_axes_locatable(ax)

            histX = divider.append_axes("top", 1.2, pad=0.2, sharex=ax)
            histY = divider.append_axes("right", 1.2, pad=0.2, sharey=ax)


            n, bins, patches = histX.hist(group['median_dsi'], bins=nbins, range=limits, normed=True, color="#000000", alpha=0.66, edgecolor='none')
            plt.setp(histX.get_xticklabels(), visible=False)
            histX.set_ylabel("Frequency")

            smoothedbins = np.arange(limits[0], limits[1], 0.01)
            kde = gaussian_kde(group['median_dsi'])
            kde.set_bandwidth(kde.factor * 0.75)
            y = kde.evaluate(smoothedbins)
            #y = y / y.max() * n.max()
            histX.plot(smoothedbins, y, 'k-')
            histX.locator_params('y', nbins=3)

            n, bins, patches = histY.hist(group['pooled_dsi'], bins=nbins, range=limits, normed=True, color="#000000", alpha=0.66, orientation='horizontal', edgecolor='none')
            plt.setp(histY.get_yticklabels(), visible=False)
            histY.set_xlabel("Frequency")
            kde = gaussian_kde(group['pooled_dsi'])
            kde.set_bandwidth(kde.factor * 0.75)
            y = kde.evaluate(smoothedbins)
            #y = y / y.max() * n.max()
            histY.plot(y, smoothedbins, 'k-')
            histY.locator_params('x', nbins=3)

            ax.set_xlim(*limits)
            ax.set_ylim(*limits)

            fig.suptitle(code, y=0.93)
            pdf.savefig(bbox_inches='tight')
            plt.close()

def dump_df(df, prefix):
    df.to_csv(prefix + '.csv', index=False, encoding='utf-8', float_format="%.10f")
    df.to_pickle(prefix + '.pkl')

def read_df_if_exists(prefix, read=True):
    if os.path.isfile(prefix + '.pkl'):
        if read:
            return pd.read_pickle(prefix + '.pkl')
        else:
            return os.path.getmtime(prefix + '.pkl')
    else:
        return None

class InvalidArgumentException(BaseException):
    pass

def main():
    import argparse
    import difflib

    parser = argparse.ArgumentParser(description='Process paired-end Illumina reads for combinatorial degron profiling')
    parser.add_argument('fastq', nargs='*', help='FASTQ files containing the reads to process. If more than two FASTQ files are given, every two consecutive files are assumed to contained paired-end reads and will be processed together. Read counts from all files will be aggregated for DSI calculation.')
    parser.add_argument('-o', '--outdir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('--fastqc', required=False, default='fastqc', help='Path to the fastq executable. If not given, fastqc will be assumed to be in PATH')
    parser.add_argument('--bowtie', required=False, default='bowtie2', help='Path to the bowtie2 executable. If not given, bowtie2 will be assumed to be in PATH')
    parser.add_argument('--phix-index', required=True, help='Path to the bowtie index of the PhiX genome.')
    parser.add_argument('--pear', required=False, default='pear', help='Path to the PEAR binary. If not given, pear will be assumed to be in PATH')
    parser.add_argument('-t', '--threads', required=False, default=1, help='Number of threads to use',  type=int)
    parser.add_argument('-c', '--configuration', required=True, help='JSON configuration file.')
    parser.add_argument('-r', '--resume', required=False, action='store_true', help='Resume aborted run? If intermediate files are found, they will be reused instead.')
    parser.add_argument('-l', '--log-level', required=False, default='INFO', choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'])
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version, help='Print version and exit')

    parser.add_argument('--xkcd', action='store_true', required=False, help=argparse.SUPPRESS)

    args = parser.parse_args()

    starttime = time.monotonic()

    os.makedirs(args.outdir, exist_ok=True)
    logging.basicConfig(filename=os.path.join(args.outdir, 'log.txt'), filemode='w', format='%(levelname)s:%(asctime)s:%(message)s', level=getattr(logging, args.log_level))

    logging.info("%s version %s" % (parser.prog, version))
    logging.info(" ".join(sys.argv))

    if len(args.fastq) % 2:
        logging.error("Uneven number of FASTQ files given, I don't know what to do. Aborting.")
        raise InvalidArgumentException("Uneven number of FASTQ files")

    # do this right away to make the user immediately aware of any exceptions that might occur due to
    # a malformed config file
    config = json.load(open(args.configuration))
    experiments = []
    for k,v in config['experiments'].items():
        exp = PyExperiment(k, v)
        experiments.append(exp)
        logging.debug(json.dumps(exp, indent=4, cls=PyExperimentJSONEncoder))

    insert_mismatches = config.get('insert_mismatches', 0)
    if config.get('barcode_match_algo', 'hamming').lower() == 'seqlev':
        counter = PySeqlevReadCounter(experiments, insert_mismatches, config.get('barcode_mismatches', 0))
    else:
        counter = PyHammingReadCounter(experiments, insert_mismatches, config.get('barcode_length', 0), config.get('barcode_mismatches', 0))
        logging.debug("Unique forward barcodes: %s" % json.dumps(counter.unique_forward_barcodes, indent=4))
        logging.debug("Unique reverse barcodes: %s" % json.dumps(counter.unique_reverse_barcodes, indent=4))

    fqcoutdir = os.path.join(args.outdir,'fastqc')
    os.makedirs(fqcoutdir, exist_ok=True)
    intermediate_outdir = os.path.join(args.outdir, 'intermediates')
    os.makedirs(intermediate_outdir, exist_ok=True)

    for fw, rev in zip(args.fastq[::2], args.fastq[1::2]):
        logging.info("Processing FASTQ pair %s and %s" % (fw, rev))
        if not args.resume or not os.path.isfile(os.path.join(fqcoutdir, os.path.splitext(os.path.basename(fw))[0]) + "_fastqc.html") or not os.path.isfile(os.path.join(fqcoutdir, os.path.splitext(os.path.basename(rev))[0]) + "_fastqc.html"):
            args.resume = False
            if exec_with_logging([args.fastqc, '--outdir=%s' % fqcoutdir] + [fw, rev], "fastqc"):
                return 1
        else:
            logging.info("Found fastqc output and resume is requested, continuing")

        fqnames = None
        bowtiefqname = None
        mergedfqname = None
        fastq = [os.path.basename(fq) for fq in (fw, rev)]
        if len(fastq[0]) == len(fastq[1]):
            s = difflib.SequenceMatcher(a=fastq[0], b=fastq[1])
            if s.ratio() >= 2 * (len(fastq[0]) - 1) / (2 * len(fastq[1])):
                b = s.get_matching_blocks()
                if len(b) == 2:
                    b.prepend((0,0,0))
                template = '%s{0}%s' % (fastq[0][b[0][0]:b[0][0]+b[0][2]], fastq[0][b[1][0]:b[1][0]+b[1][2]])
                fqnames = [os.path.join(intermediate_outdir, template.format(i)) for i in range(1,3)]
                bowtiefqname = os.path.join(intermediate_outdir, template.format('%'))
                mergedfqname = fastq[0][b[0][0]:b[0][0]+b[0][2]]

        if not fqnames:
            fqnames = [os.path.join(intermediate_outdir, 'sequence_%d.fastq' % i) for i in range(1,3)]
            bowtiefqname = os.path.join(intermediate_outdir, 'sequence_%.fastq')
            mergedfqname = 'sequence'

        bowtieout = os.path.join(intermediate_outdir, 'phix_alignment_summary.txt')
        bowtiesam = os.path.join(intermediate_outdir, 'phix_alignment.sam')
        bowtiemetrics = os.path.join(intermediate_outdir, 'phix_alignment_metrics.txt')

        if not args.resume or not os.path.isfile(bowtieout) or not os.path.isfile(bowtiesam) or not os.path.isfile(bowtiemetrics) or not os.path.isfile(fqnames[0]) or not os.path.isfile(fqnames[1]):
            args.resume = False
            if exec_with_logging([args.bowtie, '-p', str(args.threads), '--local', '--un-conc', bowtiefqname, '-x', args.phix_index, '-1', fw, '-2', rev, '-S', bowtiesam, '--no-unal', '--met-file', bowtiemetrics], "bowtie2", err=bowtieout):
                return 1
        else:
            logging.info("Found bowtie2 output and resume is requested, continuing")

        mergedfqpath = os.path.join(intermediate_outdir, mergedfqname)
        mergedfqname += '.assembled.fastq'
        pearout = os.path.join(intermediate_outdir, 'pear_summary.txt')
        if not args.resume or not os.path.isfile(pearout) or not os.path.isfile(os.path.join(intermediate_outdir, mergedfqname)):
            args.resume = False
            if exec_with_logging([args.pear, '-j', str(args.threads), '-f', fqnames[0], '-r', fqnames[1], '-o', mergedfqpath], "PEAR", out=pearout):
                return 1
        else:
            logging.info("Found PEAR output and resume is requested, continuing")

        unmatcheddir = os.path.join(args.outdir, "%s_unmapped" % mergedfqname)
        os.makedirs(unmatcheddir, exist_ok=True)

        prefixes = {}
        rawcountsprefixes = {}
        have_counts = True
        for e in experiments:
            if not len(e.name):
                prefixes[e] = ''
            else:
                prefixes[e] = "%s_" % e.name
            rawcountsprefixes[e] = os.path.join(args.outdir, "%sraw_counts" % prefixes[e])
            ex = read_df_if_exists(rawcountsprefixes[e], False)
            if ex is None or ex < os.path.getmtime(args.configuration):
                have_counts = False

        if args.resume and have_counts:
            logging.info("Found raw counts for all experiments and resume is requested, continuing")
        else:
            args.resume = False
            with TimeLogger("counting reads"):
                counter.countReads(os.path.join(intermediate_outdir, mergedfqname), os.path.join(unmatcheddir, "unmapped"), args.threads)
            logging.info("{:,d} total reads".format(counter.read))
            logging.info("{:,d} counted reads".format(counter.counted))
            logging.info("{:,d} unmatched reads, thereof {:,d} reads that could not be matched to an insert, {:,d} reads without a barcode, {:,d} reads that could not be matched to a named insert".format(counter.unmatched_total, counter.unmatched_insert, counter.unmatched_barcodes, counter.unmatched_insert_sequence))

    for e in experiments:
        if args.resume and have_counts:
            counts = read_df_if_exists(rawcountsprefixes[e])
        else:
            args.resume = False
            counts = e.counts_df
            if len(e.sorted_cells):
                with TimeLogger("normalizing counts for experiment %s" % e.name):
                    counts = normalizeCounts(counts, e.sorted_cells_df)
            counts['translation'] = pd.Series([str(Bio.Seq.Seq(str(x.sequence), Bio.Alphabet.generic_dna).translate()) for x in counts.itertuples()])
            dump_df(counts, rawcountsprefixes[e])

        if e.dsi != DSIS.noDSI:
            dspecfile = os.path.join(intermediate_outdir, "dspec.pkl")
            byaafile = os.path.join(args.outdir, "%sDSIs_byaa" % prefixes[e])
            bynucfile = os.path.join(args.outdir, "%sDSIs_bynuc" % prefixes[e])

            dspec = False
            dsi_byaa = False
            dsi_bynuc = False
            if args.resume and os.path.isfile(dspecfile):
                dspec = pickle.load(open(dspecfile, 'r+b'))
            if args.resume:
                dsi_byaa = read_df_if_exists(byaafile)
                dsi_bynuc = read_df_if_exists(bynucfile)
            if dspec is not False and dsi_byaa is not False and dsi_bynuc is not False:
                logging.info("Found DSI data for experiment %s and resume is requested, continuing" % e.name)
            else:
                args.resume = False
                dspec = getDSISpec(e)
                if dspec:
                    with TimeLogger("calculating DSIs for experiment %s" % e.name):
                        dsi_byaa, dsi_bynuc = getDSI(counts, dspec)

                    dump_df(dsi_byaa, os.path.join(args.outdir, "%sDSIs_byaa" % prefixes[e]))
                    dump_df(dsi_bynuc, os.path.join(args.outdir, "%sDSIs_bynuc" % prefixes[e]))
                    pickle.dump(dspec, open(dspecfile, 'w+b'), 3)
            if dspec:
                if args.xkcd:
                    plt.xkcd()

                ofile = os.path.join(args.outdir, "%sDSIs_byaa_cor.pdf" % prefixes[e])
                with TimeLogger("plotting DSI correlations for experiment %s into %s" % (e.name, ofile), "finished plotting DSI correlations"):
                    plot_correlations(dsi_byaa, dspec, (1, counts[dspec.dsicol].cat.categories.size), ofile)
                ofile = os.path.join(args.outdir, "%sDSIs_byaa_cor_nostop.pdf" % prefixes[e])
                with TimeLogger("plotting DSI correlations for experiment %s into %s" % (e.name, ofile), "finished plotting DSI correlations"):
                    plot_correlations(dsi_byaa[~dsi_byaa['translation'].str.contains('*', regex=False)], dspec, (1, counts[dspec.dsicol].cat.categories.size), ofile)

                histogramdir = os.path.join(args.outdir, "readcount_histograms")
                os.makedirs(histogramdir, exist_ok=True)
                for q in np.linspace(0.9, 1, 11):
                    ofile = os.path.join(histogramdir, "%sDSIs_bynuc_readcounts_%.2f_quantile.pdf" % (prefixes[e], q))
                    with TimeLogger("plotting read count histograms for experiment %s and quantile %.2f" % (e.name, q), "finished plotting read count histograms"):
                        plot_histograms(dsi_bynuc, dspec, ofile, q)
                        plot_histograms(dsi_byaa, dspec, os.path.join(histogramdir, "%sDSIs_byaa_readcounts_%.2f_quantile.pdf" % (prefixes[e], q)), q)

                ofile = os.path.join(args.outdir, "%sbynuc_countplots.pdf" % prefixes[e])
                with TimeLogger("plotting read count profiles for experiment %s into %s" % (e.name, ofile), "finished plotting read count profiles"):
                    plot_profiles(counts, [dspec.groupby, dspec.seqcol], dspec, ofile)
                ofile = os.path.join(args.outdir, "%sbyaa_countplots.pdf" % prefixes[e])
                with TimeLogger("plotting read count profiles for experiment %s into %s" % (e.name, ofile), "finished plotting read count profiles"):
                    plot_profiles(counts.groupby([dspec.groupby, dspec.dsicol, 'translation'])['normalized_counts'].sum().reset_index(), [dspec.groupby, 'translation'], dspec, ofile)

    stoptime = time.monotonic()
    logging.info("%s finished after %s" % (parser.prog, formatTime(stoptime - starttime)))
    return 0
