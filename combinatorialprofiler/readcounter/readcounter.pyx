# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref

from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.memory cimport shared_ptr
from libc.stdint cimport *

import pandas as pd



cdef extern from "cReadCounter.h" nogil:
    ctypedef unordered_map[string, unordered_set[string]] InsertSet
    ctypedef unordered_map[string, unordered_map[string, vector[string]]] BarcodeSet

    ctypedef unordered_map[string, string] SequenceSet
    ctypedef unordered_map[string, unordered_map[string, double]] SortedCellCounts
    ctypedef unordered_map[pair[string, string], unordered_map[string, uint64_t]] Counts
    ctypedef unordered_map[string, vector[string]] UniqueBarcodes

    cpdef enum NDSIS:
        noNDSI = 0
        forward = 1
        reverse = 2

    cdef cppclass Experiment:
        Experiment() except+
        Experiment(string) except+

        string name
        string insert
        SequenceSet fwBarcodeSet
        SequenceSet revBarcodeSet
        SequenceSet namedInserts
        NDSIS ndsi
        SortedCellCounts sortedCells
        Counts counts

    UniqueBarcodes makeUnique(unordered_set[string], uint16_t)

    cdef cppclass ReadCounter:
        ReadCounter(vector[Experiment*], uint16_t mismatches, uint16_t unique_barcode_length, float allowed_barcode_mismatches) except+
        void countReads(const string&, const string&, int)
        uint16_t minimumUniqueBarcodeLength()
        uint16_t allowedMismatches()
        float allowedBarcodeMismatches()
        uint64_t read()
        uint64_t unmatchedTotal()
        uint64_t unmatchedInsert()
        uint64_t unmatchedBarcodes()
        uint64_t unmatchedInsertSequence()
        uint64_t counted()
        uint64_t written()

        unordered_map[string, UniqueBarcodes] uniqueForwardBarcodes()
        unordered_map[string, UniqueBarcodes] uniqueReverseBarcodes()

cdef class PyExperiment:
    cdef Experiment *_exprmnt

    def __cinit__(self, str name, dict d=None):
        self._exprmnt = new Experiment(name.encode())

    def __dealloc__(self):
        del self._exprmnt

    def __init__(self, str name, dict d=None):
        if d:
            self.fromDict(d)

    def fromDict(self, d):
        if 'insert' in d:
            self._exprmnt.insert = d['insert'].encode()
        else:
            raise RuntimeError("Insert sequence missing from experiment %s" % self._exprmnt.name)
        if 'barcodes_fw' in d:
            for k, v in d['barcodes_fw'].items():
                self._exprmnt.fwBarcodeSet[v.encode()] = k.encode()
        if 'barcodes_rev' in d:
            for k,v in d['barcodes_rev'].items():
                self._exprmnt.revBarcodeSet[v.encode()] = k.encode()
        if 'named_inserts' in d:
            for k,v in d['named_inserts'].items():
                self._exprmnt.namedInserts[v.encode()] = k.encode()
        if 'ndsi' in d:
            haveNdsi = False
            if d['ndsi'] == 'forward':
                self._exprmnt.ndsi = forward
                haveNdsi = True
            elif d['ndsi'] == 'reverse':
                self._exprmnt.ndsi = reverse
                haveNdsi = True
            else:
                self._exprmnt.ndsi = noNDSI
            if haveNdsi and 'sortedcells' in d:
                for fw, revs in d['sortedcells'].items():
                    for r, v in revs.items():
                        self._exprmnt.sortedCells[fw.encode()][r.encode()] = v
    def toDict(self):
        d = {}
        d['insert'] = self._exprmnt.insert
        d['barcodes_fw'] = self._exprmnt.fwBarcodeSet
        d['barcodes_rev'] = self._exprmnt.revBarcodeSet
        d['named_inserts'] = self._exprmnt.namedInserts
        if self._exprmnt.ndsi == forward:
            d['ndsi'] = 'forward'
        elif self._exprmnt.ndsi == reverse:
            d['ndsi'] = 'reverse'
        d['sortedcells'] = self._exprmnt.sortedCells
        return d

    def __repr__(self):
        return "PyExperiment %s:\n%s" % (self._exprmnt.name, self.toDict())

    @property
    def counts_df(self):
        countsdict = self.counts
        insertsdict = self.named_inserts
        experiment = []
        barcode_fw = []
        barcode_rev = []
        sequence = []
        counts = []

        for codes, seqs in countsdict.items():
            for seq, cnts in seqs.items():
                experiment.append(self.name)
                barcode_fw.append(codes[0])
                barcode_rev.append(codes[1])
                sequence.append(seq)
                counts.append(cnts)
        df = pd.DataFrame()
        df['experiment'] = pd.Series(experiment, dtype='category')
        df['barcode_fw'] = pd.Categorical(barcode_fw, categories=self.forward_barcodes.values())
        df['barcode_rev'] = pd.Categorical(barcode_rev, categories=self.reverse_barcodes.values())
        df['sequence'] = sequence
        df['counts'] = counts

        if len(insertsdict):
            df['named_insert'] = [insertsdict[s] for s in sequence]

        return df

    @property
    def sorted_cells_df(self):
        countsdict = self.sorted_cells
        barcode_fw = []
        barcode_rev = []
        counts = []
        for fw, cc in countsdict.items():
            for rev, c in cc.items():
                barcode_fw.append(fw)
                barcode_rev.append(rev)
                counts.append(c)
        return pd.Series(counts, index=pd.MultiIndex.from_arrays((barcode_fw, barcode_rev), names=('barcode_fw', 'barcode_rev')))

    @property
    def name(self):
        return self._exprmnt.name

    @property
    def insert(self):
        return self._exprmnt.insert

    @property
    def forward_barcodes(self):
        return self._exprmnt.fwBarcodeSet

    @property
    def reverse_barcodes(self):
        return self._exprmnt.revBarcodeSet

    @property
    def named_inserts(self):
        return self._exprmnt.namedInserts

    @property
    def ndsi(self):
        return NDSIS(self._exprmnt.ndsi)

    @property
    def sorted_cells(self):
        return self._exprmnt.sortedCells

    @property
    def counts(self):
        return self._exprmnt.counts

def make_unique(barcodes, int minlength=0):
    cdef unordered_set[string] codes;
    for c in barcodes:
        codes.insert(c.encode())
    return makeUnique(codes, minlength);

cdef class PyReadCounter:
    cdef ReadCounter *_rdcntr
    cdef list _exprmnts

    def __cinit__(self, list experiments, int mismatches=1, int minlength=0, float barcode_mismatches=0):
        cdef vector[Experiment*] vec
        cdef PyExperiment exp;
        self._exprmnts = experiments
        for exp in self._exprmnts:
            vec.push_back((<PyExperiment>exp)._exprmnt)
        self._rdcntr = new ReadCounter(vec, mismatches, minlength, barcode_mismatches)

    def __dealloc__(self):
        del self._rdcntr

    def countReads(self, unicode fpath, unicode unmatchedpattern, threads=1):
        self._rdcntr.countReads(fpath.encode(), unmatchedpattern.encode(), threads)
        assert self.read == self.counted + self.unmatched_total
        assert self.unmatched_total == self.unmatched_insert + self.unmatched_barcodes + self.unmatched_insert_sequence
        assert self.unmatched_total == self.written

    @property
    def allowed_mismatches(self):
        return self._rdcntr.allowedMismatches()

    @property
    def allowed_barcode_mismatches(self):
        return self._rdcntr.allowedBarcodeMismatches()

    @property
    def read(self):
        return self._rdcntr.read()

    @property
    def counted(self):
        return self._rdcntr.counted()

    @property
    def unmatched_total(self):
        return self._rdcntr.unmatchedTotal()

    @property
    def unmatched_insert(self):
        return self._rdcntr.unmatchedInsert()

    @property
    def unmatched_barcodes(self):
        return self._rdcntr.unmatchedBarcodes()

    @property
    def unmatched_insert_sequence(self):
        return self._rdcntr.unmatchedInsertSequence()

    @property
    def written(self):
        return self._rdcntr.written()

    @property
    def unique_forward_barcodes(self):
        return self._rdcntr.uniqueForwardBarcodes()

    @property
    def unique_reverse_barcodes(self):
        return self._rdcntr.uniqueReverseBarcodes()

