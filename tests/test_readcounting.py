import random
import time

from combinatorialprofiler.readcounter import PyExperiment, PyReadCounter, make_unique

class FastQCreator:
    nucleotides = ['A', 'T', 'C', 'G']
    ttable = str.maketrans('ATGC', 'TACG')

    @classmethod
    def getnuc(cls, length):
        return "".join([random.choice(cls.nucleotides) for n in range(length)])

    @classmethod
    def revcompl(cls, seq):
        return seq.translate(cls.ttable)[::-1]

    def __init__(self, tmpdir, ninserts=2):
        self.totalreads = random.randint(10000, 20000)
        self.fwcodes = set()
        self.revcodes = set()
        self.inserts = set()
        self.insertlength = random.randint(100, 300)
        self.varlength = random.randint(6, 10)

        self.varcounts = [{} for i in range(ninserts)]
        self.vartotalcounts = [0] * ninserts
        self.unique_inserts = set()

        self.file = tmpdir.join("%i.fastq" % time.time())

        ncodes = random.randint(4, 10)
        while len(self.fwcodes) != ncodes:
            self.fwcodes.add(self.getnuc(random.randint(6,20)))

        ncodes = random.randint(4, 10)
        while len(self.revcodes) != ncodes:
            self.revcodes.add(self.getnuc(random.randint(6,20)))
        self.fwcodes = list(self.fwcodes)
        self.revcodes = list(self.revcodes)

        while len(self.inserts) != ninserts:
            nucs = []
            for j in range(2):
                nucs.append(self.getnuc(self.insertlength - self.varlength // 2))
            self.inserts.add(("N" * self.varlength).join(nucs))
        self.inserts = list(self.inserts)

        with self.file.open('w') as fq:
            for i in range(self.totalreads):
                ins = random.randrange(0, len(self.inserts))
                varins = self.getnuc(self.varlength)
                insert = self.inserts[ins].replace("N" * self.varlength, varins)
                self.unique_inserts.add(varins)

                fwi = random.randrange(0, len(self.fwcodes))
                revi = random.randrange(0, len(self.revcodes))
                fwcode = self.fwcodes[fwi]
                revcode = self.revcodes[revi]
                self.varcounts[ins].setdefault((str(fwi), str(revi)), {})
                self.varcounts[ins][(str(fwi), str(revi))][varins] =  self.varcounts[ins][(str(fwi), str(revi))].get(varins, 0) + 1
                self.vartotalcounts[ins] += 1
                print("@%i    %s   %i    %s   %s" % (i, fwcode, ins, varins, revcode), file=fq)
                print("".join([fwcode, insert, self.revcompl(revcode)]), file=fq)
                print("+", file=fq)
                print("I" * (len(fwcode) + len(insert) + len(revcode)), file=fq)

class ReducedCountsHolder:
    pass

def reduce_barcodes(fq, cindex, tokeep_fw, tokeep_rev):
    rc = ReducedCountsHolder()
    rc.totalcounts = [fq.vartotalcounts[cindex], 0, 0, 0]
    rc.counts = [fq.varcounts[cindex], {}, {}, {}]

    if isinstance(tokeep_fw, float) or isinstance(tokeep_fw, int):
        rc.fwcodes_tokeep = {str(i) : fq.fwcodes[i] for i in random.sample(range(len(fq.fwcodes)), int(tokeep_fw * len(fq.fwcodes)))}
    else:
        rc.fwcodes_tokeep = tokeep_fw
    rc.fwcodes_toremove = {str(i) : fq.fwcodes[i] for i in range(len(fq.fwcodes)) if str(i) not in rc.fwcodes_tokeep}
    if isinstance(tokeep_rev, float) or isinstance(tokeep_rev, int):
        rc.revcodes_tokeep = {str(i) : fq.revcodes[i] for i in random.sample(range(len(fq.revcodes)), int(tokeep_rev * len(fq.revcodes)))}
    else:
        rc.revcodes_tokeep = tokeep_rev
    rc.revcodes_toremove = {str(i): fq.revcodes[i] for i in range(len(fq.revcodes)) if str(i) not in rc.revcodes_tokeep}

    unique_forward = make_unique(rc.fwcodes_tokeep.values(), 1)
    unique_reverse = make_unique(rc.revcodes_tokeep.values(), 1)

    for codes in list(rc.counts[-4].keys()):
        fwcode = None
        revcode = None
        if codes[0] in rc.fwcodes_toremove:
            for c, u in unique_forward.items():
                for uc in u:
                    if rc.fwcodes_toremove[codes[0]].startswith(uc):
                        fwcode = c
                        break
                if fwcode:
                    for j, c in rc.fwcodes_tokeep.items():
                        if c == fwcode:
                            fwcode = j
                            break
                    break
        else:
            fwcode = codes[0]
        if codes[1] in rc.revcodes_toremove:
            for c, u in unique_reverse.items():
                for uc in u:
                    if rc.revcodes_toremove[codes[1]].startswith(uc):
                        revcode = c
                        break
                if revcode:
                    for j, c in rc.revcodes_tokeep.items():
                        if c == revcode:
                            revcode = j
                            break
                    break
        else:
            revcode = codes[1]

        index = -4
        ncodes = (fwcode, revcode)
        if fwcode not in rc.fwcodes_tokeep:
            ncodes = ('', ncodes[1])
            index += 1
        if revcode not in rc.revcodes_tokeep:
            ncodes = (ncodes[0], '')
            index += 2
        if ncodes not in rc.counts[-4]:
            ndict = rc.counts[index].setdefault(ncodes, {})
            for i, c in rc.counts[-4][codes].items():
                ndict[i] = ndict.get(i, 0) + c
                rc.totalcounts[index] += c
                rc.totalcounts[-4] -= c
            del rc.counts[-4][codes]
        elif ncodes[0] != codes[0] or ncodes[1] != codes[1]:
            for i, c in rc.counts[-4][codes].items():
                rc.counts[-4][ncodes][i] = rc.counts[-4][ncodes].get(i,0) + c
            del rc.counts[-4][codes]
    return rc

def make_barcodes_dict(codes):
    return {str(i): c for i,c in enumerate(codes)}

def simplecounts(tmpdir, mismatches):
    fq = FastQCreator(tmpdir)
    exps = []
    for i in enumerate(fq.inserts):
        d = {}
        d['insert'] = i[1]
        d['barcodes_fw'] = make_barcodes_dict(fq.fwcodes)
        d['barcodes_rev'] = make_barcodes_dict(fq.revcodes)
        exps.append(PyExperiment(str(i[0]), d))
    counter = PyReadCounter(exps, mismatches, 1)

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)

    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth
    assert counter.allowed_mismatches == mismatches

def test_simplecounts(tmpdir):
    random.seed(42)
    simplecounts(tmpdir, 1)

def test_simplecounts_nomismatches(tmpdir):
    simplecounts(tmpdir, 0)

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
    d = {'insert': fq.inserts[0], 'barcodes_fw': make_barcodes_dict(fq.fwcodes), 'barcodes_rev': make_barcodes_dict(fq.revcodes), 'named_inserts': nins}
    exp = PyExperiment('test', d)
    counter = PyReadCounter([exp], minlength=1)

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
    assert counter.allowed_mismatches == 1

def test_no_reverse_codes(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=2)
    rc = reduce_barcodes(fq, 1, 1, 0.7)
    counts = [fq.varcounts[0]] + rc.counts[::2]
    fwcodes = make_barcodes_dict(fq.fwcodes)
    exps = [
        PyExperiment("fullset", {'insert': fq.inserts[0], 'barcodes_fw':  fwcodes, 'barcodes_rev': make_barcodes_dict(fq.revcodes)}),
        PyExperiment("reduced", {'insert': fq.inserts[1], 'barcodes_fw': fwcodes, 'barcodes_rev': rc.revcodes_tokeep}),
        PyExperiment("norevcodes", {'insert': fq.inserts[1], 'barcodes_fw': fwcodes})
    ]
    counter = PyReadCounter(exps, minlength=1)
    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, counts):
        assert val.counts == truth, "Counts don't match with expected values for experiment %s" % val.name
    assert counter.allowed_mismatches == 1

def test_no_forward_codes(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=2)
    rc = reduce_barcodes(fq, 1, 0.7, 1)
    counts = [fq.varcounts[0]] + rc.counts
    revcodes = make_barcodes_dict(fq.revcodes)
    exps = [
        PyExperiment("fullset", {'insert': fq.inserts[0], 'barcodes_rev':  revcodes, 'barcodes_fw': make_barcodes_dict(fq.fwcodes)}),
        PyExperiment("reduced", {'insert': fq.inserts[1], 'barcodes_rev': revcodes, 'barcodes_fw': rc.fwcodes_tokeep}),
        PyExperiment("nofwcodes", {'insert': fq.inserts[1], 'barcodes_rev': revcodes})
    ]
    counter = PyReadCounter(exps, minlength=1)
    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, counts):
        assert val.counts == truth, "Counts don't match with expected values for experiment %s" % val.name
    assert counter.allowed_mismatches == 1

def test_no_forward_and_reverse_codes(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=2)
    rc = reduce_barcodes(fq, 1, 0.7, 0.7)
    counts = [fq.varcounts[0]] + rc.counts
    exps = [
        PyExperiment("fullset", {'insert': fq.inserts[0], 'barcodes_fw':  make_barcodes_dict(fq.fwcodes), 'barcodes_rev': make_barcodes_dict(fq.revcodes)}),
        PyExperiment("reduced", {'insert': fq.inserts[1], 'barcodes_fw': rc.fwcodes_tokeep, 'barcodes_rev': rc.revcodes_tokeep}),
        PyExperiment("nofwcodes", {'insert': fq.inserts[1], 'barcodes_rev': rc.revcodes_tokeep}),
        PyExperiment("norevcodes", {'insert': fq.inserts[1], 'barcodes_fw': rc.fwcodes_tokeep}),
        PyExperiment("nofwrevcodes", {'insert': fq.inserts[1]})
    ]
    counter = PyReadCounter(exps)

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, counts):
        assert val.counts == truth, "Counts don't match with expected values for experiment %s" % val.name
    assert counter.allowed_mismatches == 1

def test_unmatchable_inserts(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=3)
    nunmatched = 0
    for c,i in fq.varcounts[-1].items():
        for ins, c in i.items():
            nunmatched += c

    exps = [PyExperiment(str(i), {'insert': iseq, 'barcodes_fw': make_barcodes_dict(fq.fwcodes), 'barcodes_rev': make_barcodes_dict(fq.revcodes)}) for i, iseq in enumerate(fq.inserts[:-1])]
    counter = PyReadCounter(exps, minlength=1)

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    assert counter.unmatched_insert == nunmatched
    assert counter.unmatched_insert == counter.unmatched_total
    assert counter.counted == counter.read - nunmatched

    for val, truth in zip(exps, fq.varcounts[:-1]):
        assert val.counts == truth
    assert counter.allowed_mismatches == 1

def test_unmatchable_barcodes(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=2)
    nunmatched = 0
    rc = reduce_barcodes(fq, 0, 0.8, 0.8)
    nunmatched += sum(rc.totalcounts[1:])
    fwcodes = rc.fwcodes_tokeep
    revcodes = rc.revcodes_tokeep
    for i in range(1, len(fq.inserts)):
        rc = reduce_barcodes(fq, i, fwcodes, revcodes)
        nunmatched += sum(rc.totalcounts[1:])
    exps = [PyExperiment(str(i), {'insert': iseq, 'barcodes_fw': fwcodes, 'barcodes_rev': revcodes}) for i, iseq in enumerate(fq.inserts)]
    counter = PyReadCounter(exps, minlength=1)

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    #assert counter.unmatched_barcodes == nunmatched
    #assert counter.counted == counter.read - nunmatched

    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth
    assert counter.allowed_mismatches == 1
