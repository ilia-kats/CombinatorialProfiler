import random
import string

import numpy as np, numpy.random
import pandas as pd

import Bio.Seq
import Bio.Alphabet

from combinatorialprofiler.CombinatorialProfiler import normalizeCounts, getDSI, getDSISpec
from combinatorialprofiler.readcounter import DSIS

class DFCreator:
    nucleotides = ['A', 'T', 'C', 'G']

    @classmethod
    def getName(cls, length):
        return "".join([random.choice(string.ascii_letters) for i in range(length)])

    def __init__(self):
        self.seqlen = random.randrange(6, 15, 3)
        self.seqs = set(["".join([random.choice(self.nucleotides) for n in range(self.seqlen)]) for s in range(random.randint(1000, 10000))])
        self.aaseqs = set([str(Bio.Seq.Seq(s, Bio.Alphabet.generic_dna).translate()) for s in self.seqs])

        self.dsi_aa_median = {}
        self.dsi_aa_pooled = {}
        self.dsi_nuc = {}

        normalizedreads_nuc = {}
        normalizedreads_aa = {}
        totalreads = {}
        rawreads_nuc = {}
        rawreads_aa = {}

        fwcodes = set(["".join([random.choice(string.ascii_letters) for n in range(random.randint(5, 15))]) for s in range(random.randint(4, 10))])
        revcodes = [string.ascii_uppercase[i] for i in range(random.randint(3,6))]

        sortedcells = {fwc : {code: fraction for code,fraction in zip(revcodes, np.random.dirichlet(np.ones(len(revcodes)), size=1).squeeze())} for fwc in fwcodes}

        for fw in fwcodes:
            for rev in revcodes:
                rawreads_nuc[(fw, rev)] = {}
                rawreads_aa[(fw, rev)] = {}
                normalizedreads_nuc[(fw, rev)] = {}
                normalizedreads_aa[(fw, rev)] = {}
                for seq in self.seqs:
                    if random.random() > 0.6:
                        creads = random.randint(0, 1000)
                        aa = str(Bio.Seq.Seq(seq, Bio.Alphabet.generic_dna).translate())
                        rawreads_nuc[(fw, rev)][seq] = creads
                        rawreads_aa[(fw, rev)][aa] = rawreads_aa[(fw, rev)].get(aa, 0) + creads
                        totalreads[(fw, rev)] = totalreads.get((fw, rev), 0) + creads

        for cc, seqs in rawreads_nuc.items():
            normalizedreads_nuc[cc] = {}
            for s,c in seqs.items():
                normalizedreads_nuc[cc][s] = c / totalreads[cc] * sortedcells[cc[0]][cc[1]]

        for cc, seqs in rawreads_aa.items():
            normalizedreads_aa[cc] = {}
            for s,c in seqs.items():
                normalizedreads_aa[cc][s] = c / totalreads[cc] * sortedcells[cc[0]][cc[1]]

        for fw in fwcodes:
            for s in self.seqs:
                numerator = 0
                denominator = 0
                for r, rev in enumerate(revcodes):
                    if (fw, rev) in normalizedreads_nuc and s in normalizedreads_nuc[(fw, rev)]:
                        numerator += (r / (len(revcodes) - 1)) * normalizedreads_nuc[(fw, rev)][s]
                        denominator += normalizedreads_nuc[(fw, rev)][s]
                if denominator > 0:
                    self.dsi_nuc.setdefault(fw, {})[s] = numerator / denominator
                    aa = str(Bio.Seq.Seq(s, Bio.Alphabet.generic_dna).translate())
                    self.dsi_aa_median.setdefault(fw, {}).setdefault(aa, []).append(numerator / denominator)
            for s in self.aaseqs:
                numerator = 0
                denominator = 0
                for r, rev in enumerate(revcodes):
                    if (fw, rev) in normalizedreads_aa and s in normalizedreads_aa[(fw, rev)]:
                        numerator += (r / (len(revcodes) - 1)) * normalizedreads_aa[(fw, rev)][s]
                        denominator += normalizedreads_aa[(fw, rev)][s]
                if denominator > 0:
                    self.dsi_aa_pooled.setdefault(fw, {})[s] = numerator / denominator
        for fw, nn in self.dsi_aa_median.items():
            for aa, l in nn.items():
                nn[aa] = np.median(l)

        barcode_fw = []
        barcode_rev = []
        sequence = []
        counts = []
        for cc, seqs in rawreads_nuc.items():
            for s, c in seqs.items():
                barcode_fw.append(cc[0])
                barcode_rev.append(cc[1])
                sequence.append(s)
                counts.append(c)
        self.df = pd.DataFrame()
        self.df['experiment'] = pd.Series([''] * len(barcode_fw), dtype='category')
        self.df['barcode_fw'] = pd.Categorical(barcode_fw, categories=fwcodes)
        self.df['barcode_rev'] = pd.Categorical(barcode_rev, categories=revcodes)
        self.df['sequence'] = sequence
        self.df['counts'] = counts

        barcode_fw = []
        barcode_rev = []
        counts = []
        for fw, cc in sortedcells.items():
            for rev, c in cc.items():
                barcode_fw.append(fw)
                barcode_rev.append(rev)
                counts.append(c)
        self.sortedcells_df = pd.Series(counts, index=pd.MultiIndex.from_arrays((barcode_fw, barcode_rev), names=('barcode_fw', 'barcode_rev')), name='sortedcells')

def test_dsicalc():
    random.seed(42)
    simu = DFCreator()
    ncounts = normalizeCounts(simu.df, simu.sortedcells_df)
    ncounts['translation'] = pd.Series([str(Bio.Seq.Seq(str(x.sequence), Bio.Alphabet.generic_dna).translate()) for x in ncounts.itertuples()])
    dsi_byaa, dsi_bynuc = getDSI(ncounts, getDSISpec(DSIS.reverse))
    for dsi in dsi_byaa.itertuples():
        assert round(dsi.median_dsi, 5) == round(simu.dsi_aa_median[dsi.barcode_fw][dsi.translation], 5)
        assert round(dsi.pooled_dsi, 5) == round(simu.dsi_aa_pooled[dsi.barcode_fw][dsi.translation], 5)
    for dsi in dsi_bynuc.itertuples():
        assert round(dsi.dsi, 5) == round(simu.dsi_nuc[dsi.barcode_fw][dsi.sequence], 5)
