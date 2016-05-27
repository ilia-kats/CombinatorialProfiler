#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import os.path
import subprocess
import csv

class BarcodeSet:
    ttable = str.maketrans('ATGC', 'TACG')

    def __init__(self, fpath):
        with open(fpath) as bcodefile:
            sample = bcodefile.read(1024)
            bcodefile.seek(0)
            s = csv.Sniffer()
            dialect = s.sniff(sample)
            has_header = s.has_header(sample)
            reader = csv.reader(bcodefile, dialect)
            if has_header:
                next(reader)
            codes = {}
            for r in reader:
                codes[r[0]] = r[1].upper()
            self._generate_unique(codes)
            self._make_reverse()

    def _generate_unique(self, codes):
        self.fw = {}
        for k, b in codes.items():
            bcodes = []
            for s in range(0, len(b)):
                found = False
                for key, barcode in codes.items():
                    if not key == k and barcode.find(b[s:]) != -1:
                        found = True
                        break
                if not found:
                    bcodes.append(b[s:])
            self.fw[k] = tuple(bcodes)

    def _make_reverse(self):
        self.rev = {}
        for k, v in self.fw.items():
            self.rev[k] = tuple(BarcodeSet.reverse_compl(bcode) for bcode in v)

    def reverse_compl(sequence):
        sequence.translate(BarcodeSet.ttable)[::-1]

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
    parser.add_argument('-f', '--forward-barcodes', required=False, help='File containing forward barcodes. Must be csv-like with column 1 containing the name and column 2 the barcode')
    parser.add_argument('-r', '--reverse-barcodes', required=False, help='File containing reverse barcodes. Must be csv-like with column 1 containing the name and column 2 the barcode')

    args = parser.parse_args()

    fqcoutdir = os.path.join(args.outdir,'fastqc')
    if not os.path.isdir(fqcoutdir):
        os.makedirs(fqcoutdir)
    subprocess.run([args.fastqc, '--outdir=%s' % fqcoutdir, *args.fastq])

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

    with open(os.path.join(args.outdir, 'phix_alignment_summary.txt'), 'w') as phix_summary:
        subprocess.run([args.bowtie, '-p', str(args.threads), '--local', '--un-conc', os.path.join(args.outdir, bowtiefqname), '-x', args.phix_index, '-1', args.fastq[0], '-2', args.fastq[1], '-S', os.path.join(args.outdir, 'phix_alignment.sam'), '--met-file', os.path.join(args.outdir, 'phix_alignment_metrics.txt')], stderr=phix_summary)

    with open(os.path.join(args.outdir, 'pear_summary.txt'), 'w') as pear_summary:
        subprocess.run([args.pear, '-j', str(args.threads), '-f', os.path.join(args.outdir, fqnames[0]), '-r', os.path.join(args.outdir, fqnames[1]), '-o', os.path.join(args.outdir, mergedfqname)], stdout=pear_summary)
    mergedfqname += '.assembled.fastq'

    if args.forward_barcodes:
        fwcodes = BarcodeSet(args.forward_barcodes)
    if args.reverse_barcodes:
        revcodes = BarcodeSet(args.reverse_barcodes)
