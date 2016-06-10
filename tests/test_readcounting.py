import random
import time

from combinatorialprofiler.readcounter import PyExperiment, PyReadCounter

random.seed(42)

class FastQCreator:
    nucleotides = ['A', 'T', 'C', 'G']
    ttable = str.maketrans('ATGC', 'TACG')

    @classmethod
    def getnuc(cls):
        return cls.nucleotides[random.randrange(0,4)]

    @classmethod
    def revcompl(cls, seq):
        return seq.translate(cls.ttable)[::-1]

    def __init__(self, tmpdir, ninserts=2):
        self.totalreads = random.randrange(1000, 10001)
        self.fwcodes = set()
        self.revcodes = set()
        self.inserts = set()
        self.insertlength = random.randrange(100, 301)
        self.varlength = random.randrange(6, 10)

        self.varcounts = [{} for i in range(ninserts)]
        self.unique_inserts = set()

        self.file = tmpdir.join("%i.fastq" % time.time())

        ncodes = random.randrange(1, 10)
        while len(self.fwcodes) != ncodes:
            self.fwcodes.add("".join([self.getnuc() for n in range(random.randrange(6,20))]))

        ncodes = random.randrange(1, 10)
        while len(self.revcodes) != ncodes:
            self.revcodes.add("".join([self.getnuc() for n in range(random.randrange(6,20))]))
        self.fwcodes = list(self.fwcodes)
        self.revcodes = list(self.revcodes)

        while len(self.inserts) != ninserts:
            nucs = []
            for j in range(2):
                nucs.append("".join([self.getnuc() for n in range(self.insertlength - self.varlength // 2)]))
            self.inserts.add(("N" * self.varlength).join(nucs))
        self.inserts = list(self.inserts)

        with self.file.open('w') as fq:
            for i in range(self.totalreads):
                ins = random.randrange(0, len(self.inserts))
                varins = "".join([self.getnuc() for n in range(self.varlength)])
                insert = self.inserts[ins].replace("N" * self.varlength, varins)
                self.unique_inserts.add(varins)

                fwi = random.randrange(0, len(self.fwcodes))
                revi = random.randrange(0, len(self.revcodes))
                fwcode = self.fwcodes[fwi]
                revcode = self.revcodes[revi]
                self.varcounts[ins].setdefault((str(fwi), str(revi)), {})
                self.varcounts[ins][(str(fwi), str(revi))][varins] =  self.varcounts[ins][(str(fwi), str(revi))].get(varins, 0) + 1
                print("@%i    %s   %i    %s   %s" % (i, fwcode, ins, varins, revcode), file=fq)
                print("".join([fwcode, insert, self.revcompl(revcode)]), file=fq)
                print("+", file=fq)
                print("I" * (len(fwcode) + len(insert) + len(revcode)), file=fq)

def test_simplecounts(tmpdir):
    fq = FastQCreator(tmpdir)
    exps = []
    for i in enumerate(fq.inserts):
        d = {}
        d['insert'] = i[1]
        d['barcodes_fw'] = {str(i): c for i,c in enumerate(fq.fwcodes)}
        d['barcodes_rev'] = {str(i): c for i,c in enumerate(fq.revcodes)}
        exps.append(PyExperiment(str(i[0]), d))
    counter = PyReadCounter(exps)

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)

    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth

def test_namedinserts(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=1)
    ncounts = {}
    totalcounts = 0
    unmatched = 0
    nins = {i: i for i in random.sample(fq.unique_inserts, int(0.8 * len(fq.unique_inserts)))}
    for codes, vals in fq.varcounts[0].items():
        ncounts[codes] = {}
        for ins, counts in vals.items():
            if ins in nins:
                ncounts[codes][ins] = counts
                totalcounts += counts
            else:
                unmatched += counts
    d = {'insert': fq.inserts[0], 'barcodes_fw': dict((str(i), c) for i,c in enumerate(fq.fwcodes)), 'barcodes_rev': dict((str(i), c) for i,c in enumerate(fq.revcodes)), 'named_inserts': nins}
    exp = PyExperiment('test', d)
    counter = PyReadCounter([exp])

    unmatcheddir = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatcheddir), 4)

    assert counter.read == fq.totalreads
    for codes, vals in exp.counts.items():
        for i, c in vals.items():
            assert i in nins

    assert counter.counted == totalcounts
    assert counter.unmatched_insert_sequence == unmatched
    assert counter.unmatched_total == counter.unmatched_insert_sequence
    assert exp.counts == ncounts
