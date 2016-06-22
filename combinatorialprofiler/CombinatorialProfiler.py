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

import Bio.Seq
import Bio.Alphabet

from .readcounter import PyReadCounter, PyExperiment, NDSIS
from . import version

def normalizeCounts(df, sortedcells):
    df = df.set_index(['experiment','barcode_fw', 'barcode_rev','sequence'])
    df['normalized_counts'] = (df['counts'] / df.groupby(level=['barcode_fw', 'barcode_rev'])['counts'].transform(sum)).unstack(['experiment','sequence']).mul(sortedcells, axis=0).stack(['experiment','sequence']).reorder_levels(df.index.names).dropna()
    df = df.reset_index()
    df['experiment'] = df['experiment'].astype("category")
    df.barcode_fw = df.barcode_fw.astype("category")
    df.barcode_rev = df.barcode_rev.astype("category")
    return df

class NDSISpec(object):
    pass

def getNDSISpec(e):
    nspec = NDSISpec()
    if isinstance(e, PyExperiment):
        ndsi = e.ndsi
        if len(e.named_inserts):
            nspec.seqcol = 'named_insert'
        else:
            nspec.seqcol = 'sequence'
    elif isinstance(e, NDSIS):
        ndsi = e
    else:
        raise TypeError("Unsupported type")

    if ndsi == NDSIS.noNDSI:
        return None
    elif ndsi == NDSIS.forward:
        nspec.groupby = 'barcode_rev'
        nspec.ndsicol = 'barcode_fw'
    elif ndsi == NDSIS.reverse:
        nspec.groupby = 'barcode_fw'
        nspec.ndsicol = 'barcode_rev'

    return nspec

def getNDSI(df, nspec):
    statsfuns = ['min', 'max', 'mean', 'median', 'std', 'sum']

    fractionvals = pd.Series(range(1, df[nspec.ndsicol].cat.categories.size + 1), index=sorted(df[nspec.ndsicol].cat.categories))

    df['normalized_counts_cells'] = df.set_index(nspec.ndsicol, append=True)['normalized_counts'].mul(fractionvals, level=nspec.ndsicol).reset_index(level=nspec.ndsicol, drop=True)

    groupbyl = ['experiment', nspec.groupby]
    if 'named_insert' in df.columns:
        groupbyl.append('named_insert')

    g = df.groupby(groupbyl + ['translation','sequence'])
    byseq = (g['normalized_counts_cells'].sum() / g['normalized_counts'].sum()).dropna()
    byseq.name = 'ndsi'
    byseq_stats = g.agg({'counts': statsfuns, 'normalized_counts': statsfuns})
    byseq_stats.columns = ['_'.join(col) for col in byseq_stats.columns.values]
    byseq_stats['nfractions'] = g.agg('size')

    byaa_median = byseq.groupby(level=groupbyl + ['translation']).median().dropna()
    byaa_median.name = 'median_ndsi'

    g = df.groupby(groupbyl + ['translation'])
    byaa_pooled = (g['normalized_counts_cells'].sum() / g['normalized_counts'].sum()).dropna()
    byaa_pooled.name = 'pooled_ndsi'

    byaa_stats = g.agg({'counts': statsfuns, 'normalized_counts': statsfuns})
    byaa_stats.columns = ['_'.join(col) for col in byaa_stats.columns.values]
    byaa_stats['nfractions'] = g[nspec.ndsicol].nunique()

    return (pd.concat((byaa_median, byaa_pooled, byaa_stats), axis=1).reset_index().dropna(subset=('median_ndsi', 'pooled_ndsi')), pd.concat((byseq, byseq_stats), axis=1).reset_index().dropna(subset=['ndsi']))

class PyExperimentJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if type(o) != PyExperiment:
            return super().default(o)
        else:
            return {o.name : o.toDict()}

def exec_with_logging(args, pname, out=None, err=None):
    logging.info("Starting %s" % pname)
    logging.debug(" ".join(args))
    if out:
        outf = open(out, 'w')
    else:
        outf = sys.stdout
    if err:
        errf = open(err, 'w')
    else:
        errf = sys.stderr
    ctime1 = time.monotonic()
    ret = subprocess.call(args, stdout=outf, stderr=errf)
    ctime2 = time.monotonic()
    infostr = "%s finished after %d seconds" % (pname, round(ctime2 - ctime1))
    if not ret:
        logging.info(infostr)
    else:
        logging.error("%s with returncode %d" % (infostr, ret))
    return ret

def plot_profiles(df, groupby, nspec, filename, experiment):
    logging.info("Plotting read count profiles for experiment %s into %s" % (experiment, filename))
    ctime1 = time.monotonic()
    labels = sorted(df[nspec.ndsicol].cat.categories)
    integer_map = dict([(val, i) for i, val in enumerate(labels)])
    with PdfPages(filename) as pdf:
        for (code, seq), group in df.groupby(groupby):
            fig = plt.figure(figsize=(5,3))

            xvals = group[nspec.ndsicol].map(integer_map)
            xorder = xvals.argsort()
            plot = plt.plot(xvals.iloc[xorder], group['normalized_counts'].iloc[xorder], 'k.-')
            plt.xlim(0, len(labels) - 1)
            plt.ylim(ymin=0)
            plt.xticks(range(len(labels)), labels)
            plt.title("%s %s" % (code, seq))
            plt.ylabel("normalized counts")
            pdf.savefig(bbox_inches='tight')
            plt.close()
    ctime2 = time.monotonic()
    logging.info("Finished plotting read count profiles after %d seconds" % round(ctime2 - ctime1))

def plot_correlations(df, nspec, limits, filename, experiment):
    logging.info("Plotting NDSI correlations for experiment %s into %s" % (experiment, filename))
    ctime1 = time.monotonic()
    with PdfPages(filename) as pdf:
        for code, group in df.groupby(nspec.groupby):
            c = group['median_ndsi'].corr(group['pooled_ndsi'], method='spearman')

            fig = plt.figure(figsize=(5,5))
            ax = fig.add_subplot(111)
            ax.scatter(group['median_ndsi'], group['pooled_ndsi'], s=100, c="#000000", alpha=1/3, marker='.', edgecolor='none')
            ax.set_xlim(*limits)
            ax.set_ylim(*limits)
            ax.set_title(code)
            ax.set_xlabel("median NDSI")
            ax.set_ylabel("pooled NDSI")
            ax.text(0.1, 0.9, "Spearman's $r = %.3g$" % c, transform=ax.transAxes)
            ax.set_aspect('equal', adjustable='box', anchor='C')
            pdf.savefig(bbox_inches='tight')
            plt.close()
    ctime2 = time.monotonic()
    logging.info("Finished plotting NDSI correlations after %d seconds" % round(ctime2 - ctime1))

def dump_df(df, prefix):
    df.to_csv(prefix + '.csv', index=False, encoding='utf-8', float_format="%.10f")
    df.to_pickle(prefix + '.pkl')

def read_df_if_exists(prefix, read=True):
    if os.path.isfile(prefix + '.pkl'):
        if read:
            return pd.read_pickle(prefix + '.pkl')
        else:
            return True
    else:
        return False

class InvalidArgumentException(BaseException):
    pass

def main():
    import argparse
    import difflib

    parser = argparse.ArgumentParser(description='Process paired-end Illumina reads for combinatorial degron profiling')
    parser.add_argument('fastq', nargs='*', help='FASTQ files containing the reads to process. If more than two FASTQ files are given, every two consecutive files are assumed to contained paired-end reads and will be processed together. Read counts from all files will be aggregated for NDSI calculation.')
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

    counter = PyReadCounter(experiments, config.get('insert_mismatches', 0), config.get('barcode_length', 0), config.get('barcode_mismatches', 0))
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
            if not read_df_if_exists(rawcountsprefixes[e], False):
                have_counts = False

        if args.resume and have_counts:
            logging.info("Found raw counts for all experiments and resume is requested, continuing")
        else:
            args.resume = False
            logging.info("Counting reads")
            ctime1 = time.monotonic()
            counter.countReads(os.path.join(intermediate_outdir, mergedfqname), os.path.join(unmatcheddir, "unmapped"), args.threads)
            ctime2 = time.monotonic()
            logging.info("Finished counting reads after %d seconds" % round(ctime2 - ctime1))
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
                logging.info("Normalizing counts for experiment %s" % e.name)
                counts = normalizeCounts(counts, e.sorted_cells_df)
            counts['translation'] = pd.Series([str(Bio.Seq.Seq(str(x.sequence), Bio.Alphabet.generic_dna).translate()) for x in counts.itertuples()])
            dump_df(counts, rawcountsprefixes[e])

        if e.ndsi != NDSIS.noNDSI:
            nspecfile = os.path.join(intermediate_outdir, "nspec.pkl")
            byaafile = os.path.join(args.outdir, "%sNDSIs_byaa" % prefixes[e])
            bynucfile = os.path.join(args.outdir, "%sNDSIs_bynuc" % prefixes[e])

            nspec = False
            ndsi_byaa = False
            ndsi_bynuc = False
            if args.resume and os.path.isfile(nspecfile):
                nspec = pickle.load(open(nspecfile, 'r+b'))
            if args.resume:
                ndsi_byaa = read_df_if_exists(byaafile)
                ndsi_bynuc = read_df_if_exists(bynucfile)
            if nspec is not False and ndsi_byaa is not False and ndsi_bynuc is not False:
                logging.info("Found NDSI data for experiment %s and resume is requested, continuing" % e.name)
            else:
                args.resume = False
                nspec = getNDSISpec(e)
                if nspec:
                    logging.info("Calculating NDSIs for experiment %s" % e.name)
                    ctime1 = time.monotonic()
                    ndsi_byaa, ndsi_bynuc = getNDSI(counts, nspec)
                    ctime2 = time.monotonic()
                    logging.info("Finished calculating NDSIs after %d seconds" % round(ctime2 - ctime1))

                    dump_df(ndsi_byaa, os.path.join(args.outdir, "%sNDSIs_byaa" % prefixes[e]))
                    dump_df(ndsi_bynuc, os.path.join(args.outdir, "%sNDSIs_bynuc" % prefixes[e]))
                    pickle.dump(nspec, open(nspecfile, 'w+b'), 3)
            if nspec:
                plot_profiles(counts, [nspec.groupby, nspec.seqcol], nspec, os.path.join(args.outdir, "%sbynuc_countplots.pdf" % prefixes[e]), e.name)
                plot_profiles(counts.groupby([nspec.groupby, nspec.ndsicol, 'translation'])['normalized_counts'].sum().reset_index(), [nspec.groupby, 'translation'], nspec, os.path.join(args.outdir, "%sbyaa_countplots.pdf" % prefixes[e]), e.name)

                plot_correlations(ndsi_byaa, nspec, (1, counts[nspec.ndsicol].cat.categories.size), os.path.join(args.outdir, "%sNDSIs_byaa_cor.pdf" % prefixes[e]), e.name)
                plot_correlations(ndsi_byaa[~ndsi_byaa['translation'].str.contains('*', regex=False)], nspec, (1, counts[nspec.ndsicol].cat.categories.size), os.path.join(args.outdir, "%sNDSIs_byaa_cor_nostop.pdf" % prefixes[e]), e.name)
    stoptime = time.monotonic()
    logging.info("%s finished after %.2f hours" % (parser.prog, (stoptime - starttime) / 3600))
    return 0
