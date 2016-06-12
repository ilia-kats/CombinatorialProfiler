import random
import time

from combinatorialprofiler.readcounter import PyExperiment, PyReadCounter

random.seed(42)

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
                print("@%i    %s   %i    %s   %s" % (i, fwcode, ins, varins, revcode), file=fq)
                print("".join([fwcode, insert, self.revcompl(revcode)]), file=fq)
                print("+", file=fq)
                print("I" * (len(fwcode) + len(insert) + len(revcode)), file=fq)

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
    counter = PyReadCounter(exps, mismatches)

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)

    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth
    assert counter.allowed_mismatches == mismatches

def test_simplecounts(tmpdir):
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
    assert counter.allowed_mismatches == 1

def test_no_reverse_codes(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=2)
    fq.varcounts.append({})
    revcodes_tokeep = {str(i) : fq.revcodes[i] for i in random.sample(range(len(fq.revcodes)), int(0.7 * len(fq.revcodes)))}
    revcodes_toremove = {str(i) : fq.revcodes[i] for i in range(len(fq.revcodes)) if str(i) not in revcodes_tokeep}
    fwcodes = make_barcodes_dict(fq.fwcodes)
    exps = [
        PyExperiment("fullset", {'insert': fq.inserts[0], 'barcodes_fw':  fwcodes, 'barcodes_rev': make_barcodes_dict(fq.revcodes)}),
        PyExperiment("reduced", {'insert': fq.inserts[1], 'barcodes_fw': fwcodes, 'barcodes_rev': revcodes_tokeep}),
        PyExperiment("norevcodes", {'insert': fq.inserts[1], 'barcodes_fw': fwcodes})
    ]
    counter = PyReadCounter(exps)

    for codes in list(fq.varcounts[-2].keys()):
        if codes[1] in revcodes_toremove:
            revcode = None
            for c, u in counter.unique_reverse_barcodes[fq.inserts[-1]].items():
                for uc in u:
                    if revcodes_toremove[codes[1]].startswith(uc):
                        revcode = c
                        break
                if revcode:
                    for i, c in revcodes_tokeep.items():
                        if c == revcode:
                            revcode = i
                            break
                    break
            if revcode:
                ndict = fq.varcounts[-2].setdefault((codes[0], revcode), {})
            else:
                ndict = fq.varcounts[-1].setdefault((codes[0], ''), {})
            for i, c in fq.varcounts[-2][codes].items():
                ndict[i] = ndict.get(i, 0) + c
            del fq.varcounts[-2][codes]

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth, "Counts don't match with expected values for experiment %s" % val.name
    assert counter.allowed_mismatches == 1

def test_no_forward_codes(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=2)
    fq.varcounts.append({})
    fwcodes_tokeep = {str(i) : fq.fwcodes[i] for i in random.sample(range(len(fq.fwcodes)), int(0.7 * len(fq.fwcodes)))}
    fwcodes_toremove = {str(i) : fq.fwcodes[i] for i in range(len(fq.fwcodes)) if str(i) not in fwcodes_tokeep}
    revcodes = {str(i) : fq.revcodes[i] for i in random.sample(range(len(fq.revcodes)), int(0.8 * len(fq.revcodes)))}
    revcodes = make_barcodes_dict(fq.revcodes)
    exps = [
        PyExperiment("fullset", {'insert': fq.inserts[0], 'barcodes_rev':  revcodes, 'barcodes_fw': make_barcodes_dict(fq.fwcodes)}),
        PyExperiment("reduced", {'insert': fq.inserts[1], 'barcodes_rev': revcodes, 'barcodes_fw': fwcodes_tokeep}),
        PyExperiment("nofwcodes", {'insert': fq.inserts[1], 'barcodes_rev': revcodes})
    ]
    counter = PyReadCounter(exps)

    for codes in list(fq.varcounts[-2].keys()):
        if codes[0] in fwcodes_toremove:
            fwcode = None
            for c, u in counter.unique_forward_barcodes[fq.inserts[-1]].items():
                for uc in u:
                    if fwcodes_toremove[codes[0]].startswith(uc):
                        fwcode = c
                        break
                if fwcode:
                    for i, c in fwcodes_tokeep.items():
                        if c == fwcode:
                            fwcode = i
                            break
                    break
            if fwcode:
                ndict = fq.varcounts[-2].setdefault((fwcode, codes[1]), {})
            else:
                ndict = fq.varcounts[-1].setdefault(('', codes[1]), {})
            for i, c in fq.varcounts[-2][codes].items():
                ndict[i] = ndict.get(i, 0) + c
            del fq.varcounts[-2][codes]

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth, "Counts don't match with expected values for experiment %s" % val.name
    assert counter.allowed_mismatches == 1

def test_no_forward_and_reverse_codes(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=2)
    counts = fq.varcounts
    fq.varcounts.extend([{}, {}, {}])
    fwcodes_tokeep = {str(i) : fq.fwcodes[i] for i in random.sample(range(len(fq.fwcodes)), int(0.7 * len(fq.fwcodes)))}
    fwcodes_toremove = {str(i) : fq.fwcodes[i] for i in range(len(fq.fwcodes)) if str(i) not in fwcodes_tokeep}
    revcodes = {str(i) : fq.revcodes[i] for i in random.sample(range(len(fq.revcodes)), int(0.8 * len(fq.revcodes)))}
    revcodes_tokeep = {str(i) : fq.revcodes[i] for i in random.sample(range(len(fq.revcodes)), int(0.7 * len(fq.revcodes)))}
    revcodes_toremove = {str(i): fq.revcodes[i] for i in range(len(fq.revcodes)) if str(i) not in revcodes_tokeep}
    exps = [
        PyExperiment("fullset", {'insert': fq.inserts[0], 'barcodes_fw':  make_barcodes_dict(fq.fwcodes), 'barcodes_rev': make_barcodes_dict(fq.revcodes)}),
        PyExperiment("reduced", {'insert': fq.inserts[1], 'barcodes_fw': fwcodes_tokeep, 'barcodes_rev': revcodes_tokeep}),
        PyExperiment("nofwcodes", {'insert': fq.inserts[1], 'barcodes_rev': revcodes_tokeep}),
        PyExperiment("norevcodes", {'insert': fq.inserts[1], 'barcodes_fw': fwcodes_tokeep}),
        PyExperiment("nofwrevcodes", {'insert': fq.inserts[1]})
    ]
    counter = PyReadCounter(exps)

    for codes in list(counts[-4].keys()):
        fwcode = None
        revcode = None
        if codes[0] in fwcodes_toremove:
            for c, u in counter.unique_forward_barcodes[fq.inserts[1]].items():
                for uc in u:
                    if fwcodes_toremove[codes[0]].startswith(uc):
                        fwcode = c
                        break
                if fwcode:
                    for j, c in fwcodes_tokeep.items():
                        if c == fwcode:
                            fwcode = j
                            break
                    break
        else:
            fwcode = codes[0]
        if codes[1] in revcodes_toremove:
            for c, u in counter.unique_reverse_barcodes[fq.inserts[1]].items():
                for uc in u:
                    if revcodes_toremove[codes[1]].startswith(uc):
                        revcode = c
                        break
                if revcode:
                    for j, c in revcodes_tokeep.items():
                        if c == revcode:
                            revcode = j
                            break
                    break
        else:
            revcode = codes[1]

        index = -4
        ncodes = (fwcode, revcode)
        if fwcode not in fwcodes_tokeep:
            ncodes = ('', ncodes[1])
            index += 1
        if revcode not in revcodes_tokeep:
            ncodes = (ncodes[0], '')
            index += 2
        if ncodes not in counts[-4]:
            ndict = counts[index].setdefault(ncodes, {})
            for i, c in counts[-4][codes].items():
                ndict[i] = ndict.get(i, 0) + c
            del counts[-4][codes]
        elif ncodes[0] != codes[0] or ncodes[1] != codes[1]:
            for i, c in counts[-4][codes].items():
                counts[-4][ncodes][i] = counts[-4][ncodes].get(i,0) + c
            del counts[-4][codes]

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    assert counter.counted == fq.totalreads
    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth, "Counts don't match with expected values for experiment %s" % val.name
    assert counter.allowed_mismatches == 1

def test_unmatchable_inserts(tmpdir):
    fq = FastQCreator(tmpdir, ninserts=3)
    nunmatched = 0
    for c,i in fq.varcounts[-1].items():
        for ins, c in i.items():
            nunmatched += c

    exps = [PyExperiment(str(i), {'insert': iseq, 'barcodes_fw': make_barcodes_dict(fq.fwcodes), 'barcodes_rev': make_barcodes_dict(fq.revcodes)}) for i, iseq in enumerate(fq.inserts[:-1])]
    counter = PyReadCounter(exps)

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
    fwcodes = {str(i) : fq.fwcodes[i] for i in random.sample(range(len(fq.fwcodes)), int(0.8 * len(fq.fwcodes)))}
    fwcodes_toremove = {str(i) : fq.fwcodes[i] for i in range(len(fq.fwcodes)) if str(i) not in fwcodes}
    revcodes = {str(i) : fq.revcodes[i] for i in random.sample(range(len(fq.revcodes)), int(0.8 * len(fq.revcodes)))}
    revcodes_toremove = {str(i): fq.revcodes[i] for i in range(len(fq.revcodes)) if str(i) not in revcodes}
    exps = [PyExperiment(str(i), {'insert': iseq, 'barcodes_fw': fwcodes, 'barcodes_rev': revcodes}) for i, iseq in enumerate(fq.inserts)]
    counter = PyReadCounter(exps)

    for i,cc in enumerate(fq.varcounts):
        for codes in list(cc.keys()):
            fwcode = None
            revcode = None
            if codes[0] in fwcodes_toremove:
                for c, u in counter.unique_forward_barcodes[fq.inserts[i]].items():
                    for uc in u:
                        if fwcodes_toremove[codes[0]].startswith(uc):
                            fwcode = c
                            break
                    if fwcode:
                        for j, c in fwcodes.items():
                            if c == fwcode:
                                fwcode = j
                                break
                        break
            else:
                fwcode = codes[0]
            if codes[1] in revcodes_toremove:
                for c, u in counter.unique_reverse_barcodes[fq.inserts[i]].items():
                    for uc in u:
                        if revcodes_toremove[codes[1]].startswith(uc):
                            revcode = c
                            break
                    if revcode:
                        for j, c in revcodes.items():
                            if c == revcode:
                                revcode = j
                                break
                        break
            else:
                revcode = codes[1]
            if not fwcode or not revcode:
                for j, c in cc[codes].items():
                    nunmatched += c
                del cc[codes]
            elif fwcode != codes[0] or revcode != codes[1]:
                for j, c in cc[codes].items():
                    cc[(fwcode, revcode)][j] =  cc[(fwcode, revcode)].get(j, 0) + c
                del cc[codes]

    unmatched = tmpdir.mkdir("%s_unmatched" % fq.file.basename)
    counter.countReads(str(fq.file), str(unmatched), 4)
    assert counter.read == fq.totalreads
    assert counter.unmatched_barcodes == nunmatched
    assert counter.counted == counter.read - nunmatched

    for val, truth in zip(exps, fq.varcounts):
        assert val.counts == truth
    assert counter.allowed_mismatches == 1
