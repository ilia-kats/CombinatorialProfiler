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

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('PDF')
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

def getNDSI(df, nspec):
    if nspec == NDSIS.forward:
        groupby = 'barcode_rev'
        ndsicol = 'barcode_fw'
    elif nspec == NDSIS.reverse:
        groupby = 'barcode_fw'
        ndsicol = 'barcode_rev'
    else:
        return None

    fractionvals = pd.Series(range(1, df[ndsicol].cat.categories.size + 1), index=sorted(df[ndsicol].cat.categories))

    df['normalized_counts_cells'] = df.set_index(ndsicol, append=True)['normalized_counts'].mul(fractionvals, level=1).reset_index(level=1, drop=True)

    groupbyl = ['experiment', groupby]
    if 'named_insert' in df.columns:
        groupbyl.append('named_insert')

    g = df.groupby(groupbyl + ['translation','sequence'])
    byseq = (g['normalized_counts_cells'].sum() / g['normalized_counts'].sum()).dropna()
    byseq.name = 'ndsi'

    byaa_median = byseq.groupby(level=groupbyl + ['translation']).median().dropna()
    byaa_median.name = 'median_ndsi'

    g = df.groupby(groupbyl + ['translation'])
    byaa_pooled = g['normalized_counts_cells'].sum() / g['normalized_counts'].sum().dropna()
    byaa_pooled.name = 'pooled_ndsi'
    return (groupby, ndsicol, pd.concat((byaa_median, byaa_pooled), axis=1).dropna().reset_index(), byseq.reset_index())

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
    ret = subprocess.run(args, stdout=outf, stderr=errf).returncode
    ctime2 = time.monotonic()
    infostr = "%s finished after %i seconds" % (pname, round(ctime2 - ctime1))
    if not ret:
        logging.info(infostr)
    else:
        logging.error("%s with returncode %i" % (infostr, ret))
    return ret

def main():
    import argparse
    import difflib

    parser = argparse.ArgumentParser(description='Process paired-end Illumina MiSeq reads for combinatorial degron profiling')
    parser.add_argument('fastq', nargs=2, help='Two FASTQ files containing the forward and reverse reads, respectively')
    parser.add_argument('-o', '--outdir', required=False, default=os.getcwd(), help='Output directory')
    parser.add_argument('--fastqc', required=False, default='fastqc', help='Path to the fastq executable. If not given, fastqc will be assumed to be in PATH')
    parser.add_argument('--bowtie', required=False, default='bowtie2', help='Path to the bowtie2 executable. If not given, bowtie2 will be assumed to be in PATH')
    parser.add_argument('--phix_index', required=True, help='Path to the bowtie index of the PhiX genome.')
    parser.add_argument('--pear', required=True, help='Path to the PEAR binary. If not given, pear will be assumed to be in PATH')
    parser.add_argument('-t', '--threads', required=False, default=1, help='Number of threads to use',  type=int)
    parser.add_argument('-c', '--configuration', required=True, help='JSON configuration file.')
    parser.add_argument('-r', '--resume', required=False, action='store_true', help='Resume aborted run? If intermediate files are found, they will be reused instead.')
    parser.add_argument('-l', '--log-level', required=False, default='INFO', choices=['ERROR', 'WARNING', 'INFO', 'DEBUG'])
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + version, help='Print version and exit')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    logging.basicConfig(filename=os.path.join(args.outdir, 'run_log.txt'), filemode='w', format='%(levelname)s:%(asctime)s:%(message)s', level=getattr(logging, args.log_level))

    logging.info("%s version %s" % (parser.prog, version))
    logging.info(" ".join(sys.argv))

    # do this right away to make the user immediately aware of any exceptions that might occur due to
    # a malformed config file
    config = json.load(open(args.configuration))
    experiments = []
    for k,v in config['experiments'].items():
        exp = PyExperiment(k, v)
        experiments.append(exp)
        logging.debug(json.dumps(exp, indent=4, cls=PyExperimentJSONEncoder))

    counter = PyReadCounter(experiments, config.get('insert_mismatches', 0), config.get('barcode_length', 0))

    stime = time.monotonic()

    fqcoutdir = os.path.join(args.outdir,'fastqc')
    os.makedirs(fqcoutdir, exist_ok=True)
    if not args.resume or not os.path.isdir(fqcoutdir):
        args.resume = False
        if exec_with_logging([args.fastqc, '--outdir=%s' % fqcoutdir, *args.fastq], "fastqc"):
            return 1
    else:
        logging.info("Found fastqc output and resume is requested, continuing")

    intermediate_outdir = os.path.join(args.outdir, 'intermediates')
    os.makedirs(intermediate_outdir, exist_ok=True)
    fqnames = None
    bowtiefqname = None
    mergedfqname = None
    fastq = [os.path.basename(fq) for fq in args.fastq]
    if len(fastq[0]) == len(fastq[1]):
        s = difflib.SequenceMatcher(a=fastq[0], b=fastq[1])
        if s.ratio() >= 2 * (len(fastq[0]) - 1) / (2 * len(fastq[0])):
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
        if exec_with_logging([args.bowtie, '-p', str(args.threads), '--local', '--un-conc', bowtiefqname, '-x', args.phix_index, '-1', args.fastq[0], '-2', args.fastq[1], '-S', bowtiesam, '--met-file', bowtiemetrics], "bowtie2", err=bowtieout):
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

    logging.info("Starting counting reads")
    ctime1 = time.monotonic()
    counter.countReads(os.path.join(intermediate_outdir, mergedfqname), os.path.join(unmatcheddir, "unmapped"), args.threads)
    ctime2 = time.monotonic()
    logging.info("Finished counting reads after %i seconds" % round(ctime2 - ctime1))
    logging.debug("Unique forward barcodes: %s" % json.dumps(counter.unique_forward_barcodes, indent=4))
    logging.debug("Unique reverse barcodes: %s" % json.dumps(counter.unique_reverse_barcodes, indent=4))

    for e in experiments:
        counts = e.counts_df
        if len(e.sorted_cells):
            logging.info("Normalizing counts for experiment %s" % e.name)
            counts = normalizeCounts(counts, e.sorted_cells_df)
        counts['translation'] = pd.Series([str(Bio.Seq.Seq(str(x.sequence), Bio.Alphabet.generic_dna).translate()) for x in counts.itertuples()])

        if not len(e.name):
            prefix = ''
        else:
            prefix = "%s_" % e.name
        counts.to_csv(os.path.join(args.outdir, "%sraw_counts.csv" % prefix), index=False, encoding='utf-8')

        if e.ndsi != NDSIS.noNDSI:
            logging.info("Calculating NDSIs for experiment %s" % e.name)
            ctime1 = time.monotonic()
            groupby, ndsicol, ndsi_byaa, ndsi_bynuc = getNDSI(counts, e.ndsi)
            ctime2 = time.monotonic()
            logging.info("Finished calculating NDSIs after %i seconds" % round(ctime2 - ctime1))

            ndsi_byaa.to_csv(os.path.join(args.outdir, "%sNDSIs_byaa.csv" % prefix), index=False, encoding='utf-8')
            ndsi_bynuc.to_csv(os.path.join(args.outdir, "%sNDSIs_bynuc.csv" % prefix), index=False, encoding='utf-8')
            labels = sorted(counts[ndsicol].cat.categories)
            integer_map = dict([(val, i) for i, val in enumerate(labels)])
            if 'named_insert' in counts:
                seqcol = 'named_insert'
            else:
                seqcol = 'sequence'

            logging.info("Plotting read count profiles for experiment %s" % e.name)
            ctime1 = time.monotonic()
            with PdfPages(os.path.join(args.outdir, "%scountplots.pdf" % prefix)) as pdf:
                for (code, seq), group in counts.groupby([groupby, seqcol]):
                    fig = plt.figure(figsize=(5,3))
                    plot = plt.plot(group[ndsicol].map(integer_map), group['normalized_counts'], 'ko-')
                    plt.xlim(0, len(labels) - 1)
                    plt.ylim(ymin=0)
                    plt.xticks(range(len(labels)), labels)
                    plt.title("%s %s" % (code, seq))
                    plt.ylabel("normalized counts")
                    pdf.savefig(bbox_inches='tight')
                    plt.close()
            ctime2 = time.monotonic()
            logging.info("Finished plotting read count profiles after %i seconds" % round(ctime2 - ctime1))
    return 0
