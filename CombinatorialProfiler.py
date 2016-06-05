#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import os.path
import subprocess
import csv

import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import Bio.Seq
import Bio.Alphabet

from readcounter.readcounter import PyReadCounter

def dict_merge(dicts):
    dct = dicts.pop()
    while len(dicts):
        merge_dct = dicts.pop()
        for k, v in merge_dct.items():
            if (k in dct and isinstance(dct[k], dict)
                    and isinstance(merge_dct[k], dict)):
                dict_merge(dct[k], merge_dct[k])
            else:
                dct[k] = merge_dct[k]
    return dct

def get_csv_reader(csvfile):
    sample = csvfile.read(1024)
    csvfile.seek(0)
    s = csv.Sniffer()
    dialect = s.sniff(sample)
    has_header = s.has_header(sample)
    reader = csv.reader(csvfile, dialect)
    if has_header:
        next(reader)
    return (reader, has_header)

def split_spec(spec):
    sep = spec.find(':')
    if sep <= 0:
        return ('', spec)
    else:
        return (spec[:sep], spec[sep + 1:])

def get_named_file(spec):
    if os.path.isfile(spec):
        return (None, spec)
    else:
        name, fpath = split_spec(spec)
        if not os.path.isfile(fpath):
            raise RuntimeError("Invalid file specification")
        else:
            return (name, fpath)

def get_named_reader(spec):
    name, fpath = get_named_file(spec)
    f = open(fpath, newline='')
    return (name, get_csv_reader(f)[0], f)

def generateBarcodes(codes):
    ccodes = {}
    for i, c in codes.items():
        ccodes[i] = {}
        for k,b in c.items():
            bcodes = []
            for s in range(0, len(b)):
                found = False
                for key, barcode in codes[i].items():
                    if not key == k and barcode.find(b[s:]) != -1:
                        found = True
                        break
                if not found:
                    bcodes.append(b[s:])
            ccodes[i][k] = tuple(bcodes)
    return ccodes

def readBarcodes(fpath, reverse=False):
    ttable = str.maketrans('ATGC', 'TACG')

    name, reader, f = get_named_reader(fpath)
    codes = {}
    for r in reader:
        if len(r) > 2 and not name:
            insname = r[2]
        elif name:
            insname = name
        else:
            insname = ""
        if insname not in codes:
            codes[insname] = {}
        codes[insname][r[0]] = r[1].upper()

    if reverse:
        for i, b in codes.items():
            for k,v in b.items():
                b[k] = v.translate(ttable)[::-1]

    return codes

def readCellCounts(spec, barcodes_fw, barcodes_rev):
    name, fpath = get_named_file(spec)
    f = open(fpath, newline='')
    header = get_csv_reader(f)[1]
    f.close()

    if header:
        header = 0
    else:
        header = None
    if barcodes_fw is None and barcodes_rev is None or barcodes_fw is not None and name and name not in barcodes_fw and barcodes_rev is not None and name and name not in barcodes_rev:
        return {}
    df = pd.read_csv(fpath, header=header, sep=None)
    types = df.get_dtype_counts()
    if not name or types['float64'] != df.columns.size - 1:
        def getCountsDict(cname, group, barcodes_fw, barcodes_rev, indexcols, valuecol):
            if barcodes_fw is not None and cname in barcodes_fw and set(group.iloc[:,indexcols[0]].squeeze()) == barcodes_fw[cname].keys():
                if barcodes_rev is not None and cname in barcodes_rev and set(group.iloc[:,indexcols[1]].squeeze()) != barcodes_rev[cname].keys():
                    raise RuntimeError("Barcode labels don't match sorted cells labels")
                else:
                    index = group.columns[indexcols].values.tolist()
            elif barcodes_rev is not None and cname in barcodes_rev and set(group.iloc[:,indexcols[0]].squeeze()) == barcodes_rev[cname].keys():
                index = group.columns[list(reversed(indexcols))].values.tolist()
            else:
                raise RuntimeError("Barcode labels don't match sorted cells labels")
            cgroup = group.set_index(index)
            cgroup.index.rename(['barcode_fw', 'barcode_rev'], inplace=True)
            return cgroup.loc[:,group.columns[valuecol]]

        valuecol = np.where(df.dtypes == 'float64')[0]
        if valuecol.size > 1:
            raise RuntimeError("Unrecognized sorted cells format: Multiple numeric columns")
        if valuecol[0] == 0:
            raise RuntimeError("Unrecognized sorted cells format: No barcode labels given")
        if barcodes_fw is not None and barcodes_rev is not None:
            indexcols = [valuecol[0] - 2, valuecol[0] - 1]
        else:
            df[df.columns.size] = ''
            indexcols = [valuecol - 1, df.columns.size - 1]
        if df.columns.size > valuecol[0] and not name:
            df = df.set_index(df.columns[valuecol[0] + 1])
            retdict = {}
            for cname, group in df.groupby(level=0):
                retdict[cname] = getCountsDict(cname, group, barcodes_fw, barcodes_rev, indexcols, valuecol[0])
            return retdict
        else:
            if not name:
                name = ''
            return {name: getCountsDict(name, df, barcodes_fw, barcodes_rev, indexcols, valuecol[0])}
            raise RuntimeError("No insert name given")
    else:
        df = df.set_index(df.columns[0])
        if barcodes_fw is not None and name in barcodes_fw and barcodes_rev is not None and name in barcodes_rev:
            if set(df.index) == barcodes_fw[name].keys() and set(df.columns) == barcodes_rev[name].keys():
                toreturn = df.stack()
            elif set(df.index) == barcodes_rev[name].keys() and set(df.columns) == barcodes_fw[name].keys():
                toreturn = df.T.stack()
            else:
                raise RuntimeError("Barcode labels don't match sorted cells labels")
        else:
            if min(df.shape) != 1:
                raise RuntimeError("Too many sorted cells labels for the barcodes given")
            if df.shape[1] > df.shape[0]:
                df = df.T

            if barcodes_fw is not None and name in barcodes_fw:
                tomatch = barcodes_fw
                tdict = lambda x: x.rename(columns={x.columns[0]: ''}).stack()
            elif barcodes_rev is not None and name in barcodes_rev:
                tomatch = barcodes_rev
                tdict = lambda x: x.T.rename(columns={x.columns[0]: ''}).stack()

            if set(df.index) != tomatch[name].keys():
                raise RuntimeError("Barcode labels don't match sorted cells labels")
            else:
                toreturn = tdict(df)
        toreturn.index.rename(['barcode_fw', 'barcode_rev'], inplace=True)
        return {name: toreturn}

def readNamedInserts(ins):
    inserts = {}
    name, reader, f = get_named_reader(ins)
    inserts[name] = {}
    for r in reader:
        inserts[name][r[0]] = r[1]
    return inserts

def readInserts(ins):
    inserts = {}
    if os.path.isfile(ins):
        with open(ins) as insfile:
            reader = get_csv_reader(insfile)[0]
            for r in reader:
                inserts[r[0]] = r[1]
    else:
        name, insert = split_spec(ins)
        inserts[name] = insert
    return inserts

def mergeInserts(ins, barcodes_fw, barcodes_rev):
    if len(ins) == 1 and None in ins:
        keys = []
        if barcodes_fw is not None:
            keys.extend(barcodes_fw.keys())
        if barcodes_rev is not None:
            keys.extend(barcodes_rev.keys())
        keys = set(keys)
        insseq = ins[None]
        ins = {}
        for k in keys:
            ins[k] = insseq
    return ins

def normalizeCounts(df, sortedcells):
    df = df.set_index(['insert','barcode_fw', 'barcode_rev','sequence'])
    df['normalized_counts'] = (df['counts'] / df.groupby(level=['barcode_fw', 'barcode_rev'])['counts'].transform(sum)).unstack(['insert','sequence']).mul(sortedcells, axis=0).stack(['insert','sequence']).reorder_levels(df.index.names).dropna()
    df = df.reset_index()
    df['insert'] = df['insert'].astype("category")
    df.barcode_fw = df.barcode_fw.astype("category")
    df.barcode_rev = df.barcode_rev.astype("category")
    return df

def getNDSI(df, nspec):
    def calcNDSI(group, ndsicol, fractionvals):
        calc = lambda x: sum((x * fractionvals).dropna()) / sum(x)
        g = group.set_index(ndsicol)
        pooled = calc(g.groupby(level=ndsicol)['normalized_counts'].sum().dropna())
        med = g.groupby('sequence')['normalized_counts'].aggregate(calc).median()
        return pd.Series({'median_ndsi': med, 'pooled_ndsi': pooled})

    insert, direction = split_spec(nspec)
    if insert not in df:
        raise RuntimeError("No data for NDSI calculation: %s" % insert)
    elif direction != "forward" and direction != "reverse":
        raise RuntimeError("Unrecognized direction")
    elif direction == "forward":
        groupby = 'barcode_rev'
        ndsicol = 'barcode_fw'
    else:
        groupby = 'barcode_fw'
        ndsicol = 'barcode_rev'
    fractionvals = pd.Series(range(1, df[insert][ndsicol].cat.categories.size + 1), index=sorted(df[insert][ndsicol].cat.categories))
    return (insert, groupby, ndsicol, df[insert].groupby([groupby, 'translation']).apply(calcNDSI, ndsicol, fractionvals).dropna().reset_index())

if __name__ == '__main__':
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
    parser.add_argument('-f', '--forward-barcodes', required=False, action='append', help='File containing forward barcodes. Must be csv-like with column 1 containing the name and column 2 the barcode')
    parser.add_argument('-r', '--reverse-barcodes', required=False, action='append', help='File containing reverse barcodes. Must be csv-like with column 1 containing the name and column 2 the barcode')
    parser.add_argument('-i', '--insert-sequence', required=False, action='append', help='Either an insert sequence to match against or path to a csv-like file with column 1 containing the name and column 2 the sequence. Variable region must be marked with a sequence of Ns')
    parser.add_argument('-n', '--named-inserts', required=False, action='append', help='Named sequences of allowed variable sequences within the insert. Csv-like file with the first column containing the name and the second column the sequence')
    parser.add_argument('-s', '--sorted-cells', required=False, action='append', help='Path to a CSV file containing the relative amounts of sorted cells per barcode combination')
    parser.add_argument('--ndsi', action='append', required=False, help='Whether to use forward or reverse barcodes as FACS fractions for NDSI calculation.')

    args = parser.parse_args()

    fqcoutdir = os.path.join(args.outdir,'fastqc')
    if not os.path.isdir(fqcoutdir):
        os.makedirs(fqcoutdir)
    #subprocess.run([args.fastqc, '--outdir=%s' % fqcoutdir, *args.fastq])

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
            fqnames = [template.format(i) for i in range(1,3)]
            bowtiefqname = template.format('%')
            mergedfqname = fastq[0][b[0][0]:b[0][0]+b[0][2]]

    if not fqnames:
        fqnames = ['sequence_%d.fastq' % i for i in range(1,3)]
        bowtiefqname = 'sequence_%.fastq'
        mergedfqname = 'sequence'

    #with open(os.path.join(args.outdir, 'phix_alignment_summary.txt'), 'w') as phix_summary:
        #subprocess.run([args.bowtie, '-p', str(args.threads), '--local', '--un-conc', os.path.join(args.outdir, bowtiefqname), '-x', args.phix_index, '-1', args.fastq[0], '-2', args.fastq[1], '-S', os.path.join(args.outdir, 'phix_alignment.sam'), '--met-file', os.path.join(args.outdir, 'phix_alignment_metrics.txt')], stderr=phix_summary)

    #with open(os.path.join(args.outdir, 'pear_summary.txt'), 'w') as pear_summary:
        #subprocess.run([args.pear, '-j', str(args.threads), '-f', os.path.join(args.outdir, fqnames[0]), '-r', os.path.join(args.outdir, fqnames[1]), '-o', os.path.join(args.outdir, mergedfqname)], stdout=pear_summary)
    mergedfqname += '.assembled.fastq'

    if args.forward_barcodes:
        fwcodes = generateBarcodes(dict_merge([readBarcodes(c) for c in args.forward_barcodes]))
    else:
        fwcodes = None
    if args.reverse_barcodes:
        revcodes = generateBarcodes(dict_merge([readBarcodes(c, True) for c in args.reverse_barcodes]))
    else:
        revcodes = None

    if args.named_inserts:
        ninserts = mergeInserts(dict_merge([readNamedInserts(i) for i in args.named_inserts]), fwcodes, revcodes)
    else:
        ninserts = None
    unmatcheddir = os.path.join(args.outdir, "%s_unmapped" % mergedfqname)
    if not os.path.isdir(unmatcheddir):
        os.makedirs(unmatcheddir)
    counter = PyReadCounter(dict_merge([readInserts(ins) for ins in args.insert_sequence]), fwcodes, revcodes, ninserts)
    counter.countReads(os.path.join(args.outdir, mergedfqname), os.path.join(unmatcheddir, "unmapped"), args.threads)

    df = counter.asDataFrames()

    if args.sorted_cells:
        sortedcells = dict_merge([readCellCounts(c, fwcodes, revcodes) for c in args.sorted_cells])
        for i,v in df.items():
            if i in sortedcells:
                df[i] = normalizeCounts(v, sortedcells[i])

    for i, v in df.items():
        if not len(i):
            prefix = ''
        else:
            prefix = "%s_" % i
        v['translation'] = pd.Series([str(Bio.Seq.Seq(str(x.sequence), Bio.Alphabet.generic_dna).translate()) for x in v.itertuples()])
        v.to_csv(os.path.join(args.outdir, "%sraw_counts.csv" % prefix), index=False, encoding='utf-8')

    if args.ndsi:
        for n in args.ndsi:
            ins, groupby, ndsicol, ndsi = getNDSI(df, n)
            if not len(ins):
                prefix = ''
            else:
                prefix = "%s_" % i
            ndsi.to_csv(os.path.join(args.outdir, "%sNDSIs.csv" % prefix), index=False, encoding='utf-8')
            labels = sorted(df[ins][ndsicol].cat.categories)
            integer_map = dict([(val, i) for i, val in enumerate(labels)])
            if 'named_insert' in df[ins]:
                seqcol = 'named_insert'
            else:
                seqcol = 'sequence'
            with PdfPages(os.path.join(args.outdir, "%scountplots.pdf" % prefix)) as pdf:
                for (code, seq), group in df[ins].groupby([groupby, seqcol]):
                    fig = plt.figure(figsize=(5,3))
                    plot = plt.plot(group[ndsicol].map(integer_map), group['normalized_counts'], 'ko-')
                    plt.xlim(0, len(labels) - 1)
                    plt.xticks(range(len(labels)), labels)
                    plt.title("%s %s" % (code, seq))
                    plt.ylabel("normalized counts")
                    pdf.savefig(bbox_inches='tight')
                    plt.close()

