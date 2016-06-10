import random
import time

from readcounter import PyExperiment, PyReadCounter

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
        self.totalreads = random.randrange(10000, 100001)
        self.fwcodes = set()
        self.revcodes = set()
        self.inserts = set()
        self.insertlength = random.randrange(100, 301)
        self.varlength = random.randrange(6, 10)

        self.varcounts = [{} for i in range(ninserts)]

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

def test_init(tmpdir):
    global fq
    fq = FastQCreator(tmpdir)

def test_counts(tmpdir):
    exps = []
    for i in enumerate(fq.inserts):
        d = {}
        d['insert'] = i[1]
        d['barcodes_fw'] = dict((str(i), c) for i,c in enumerate(fq.fwcodes))
        d['barcodes_rev'] = dict((str(i), c) for i,c in enumerate(fq.revcodes))
        exps.append(PyExperiment(str(i[0]), d))
    counter = PyReadCounter(exps)

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)

    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth
