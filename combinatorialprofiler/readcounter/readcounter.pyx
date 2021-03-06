# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref

from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.memory cimport shared_ptr
from libcpp cimport bool
from libc.stdint cimport *

import pandas as pd

cdef extern from "util.h" nogil:
    ctypedef unordered_map[string, vector[string]] UniqueBarcodes
    UniqueBarcodes makeUnique(unordered_set[string], uint16_t)
    size_t hamming_distance(string, string)
    size_t seqlev_distance(string, string, bool)

cdef extern from "Experiment.h" nogil:
    ctypedef unordered_map[string, unordered_set[string]] InsertSet
    ctypedef unordered_map[string, unordered_map[string, vector[string]]] BarcodeSet

    ctypedef unordered_map[string, string] SequenceSet
    ctypedef unordered_map[string, unordered_map[string, double]] SortedCellCounts
    ctypedef unordered_map[pair[string, string], unordered_map[string, uint64_t]] Counts

    cpdef enum DSIS:
        noDSI = 0
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
        DSIS dsi
        SortedCellCounts sortedCells
        Counts counts

cdef extern from "ReadCounter.h" nogil:
    cdef cppclass ReadCounter:
        ReadCounter(vector[Experiment*], uint16_t mismatches) except+
        void countReads(const string&, const string&, int)
        uint16_t allowedMismatches()
        uint64_t read()
        uint64_t unmatchedTotal()
        uint64_t unmatchedInsert()
        uint64_t unmatchedBarcodes()
        uint64_t unmatchedInsertSequence()
        uint64_t counted()
        uint64_t written()

cdef extern from "HammingReadCounter.h" nogil:
    cdef cppclass HammingReadCounter(ReadCounter):
        @staticmethod
        HammingReadCounter* getReadCounter(vector[Experiment*], uint16_t mismatches, uint16_t unique_barcode_length, float allowed_barcode_mismatches) except+
        uint16_t minimumUniqueBarcodeLength()
        float allowedBarcodeMismatches()
        unordered_map[string, UniqueBarcodes] uniqueForwardBarcodes()
        unordered_map[string, UniqueBarcodes] uniqueReverseBarcodes()

cdef extern from "SeqlevReadCounter.h" nogil:
    cdef cppclass SeqlevReadCounter(ReadCounter):
        @staticmethod
        SeqlevReadCounter* getReadCounter(vector[Experiment*], uint16_t mismatches, uint16_t barcode_mismatches) except+
        uint16_t allowedBarcodeMismatches()

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
            self._exprmnt.insert = d['insert'].upper().encode()
        else:
            raise RuntimeError("Insert sequence missing from experiment %s" % self._exprmnt.name)
        if 'barcodes_fw' in d:
            for k, v in d['barcodes_fw'].items():
                self._exprmnt.fwBarcodeSet[v.upper().encode()] = k.encode()
        if 'barcodes_rev' in d:
            for k,v in d['barcodes_rev'].items():
                self._exprmnt.revBarcodeSet[v.upper().encode()] = k.encode()
        if 'named_inserts' in d:
            for k,v in d['named_inserts'].items():
                self._exprmnt.namedInserts[v.upper().encode()] = k.encode()
        if 'dsi' in d:
            haveDsi = False
            if d['dsi'] == 'forward':
                self._exprmnt.dsi = forward
                haveDsi = True
            elif d['dsi'] == 'reverse':
                self._exprmnt.dsi = reverse
                haveDsi = True
            else:
                self._exprmnt.dsi = noDSI
            if haveDsi and 'sortedcells' in d:
                for fw, revs in d['sortedcells'].items():
                    for r, v in revs.items():
                        self._exprmnt.sortedCells[fw.encode()][r.encode()] = v
    def toDict(self):
        d = {}
        d['insert'] = self._exprmnt.insert
        d['barcodes_fw'] = self._exprmnt.fwBarcodeSet
        d['barcodes_rev'] = self._exprmnt.revBarcodeSet
        d['named_inserts'] = self._exprmnt.namedInserts
        if self._exprmnt.dsi == forward:
            d['dsi'] = 'forward'
        elif self._exprmnt.dsi == reverse:
            d['dsi'] = 'reverse'
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
        df['barcode_fw'] = pd.Categorical(barcode_fw, categories=sorted(self.forward_barcodes.values()))
        df['barcode_rev'] = pd.Categorical(barcode_rev, categories=sorted(self.reverse_barcodes.values()))
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
        return pd.Series(counts, index=pd.MultiIndex.from_arrays((barcode_fw, barcode_rev), names=('barcode_fw', 'barcode_rev')), name='sortedcells')

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
    def dsi(self):
        return DSIS(self._exprmnt.dsi)

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

def hammingDistance(str needle, str haystack):
    return hamming_distance(needle.encode(), haystack.encode())

def seqLevDistance(str needle, str haystack, at_end=False):
    return seqlev_distance(needle.encode(), haystack.encode(), at_end)

cdef class PyReadCounter:
    cdef ReadCounter *_rdcntr
    cdef list _exprmnts

    def __cinit__(self, list experiments, *args, **kwargs):
        self._exprmnts = experiments

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

cdef class PyHammingReadCounter(PyReadCounter):
    def __cinit__(self, list experiments, int mismatches=1, int minlength=0, float barcode_mismatches=0):
        cdef vector[Experiment*] vec
        cdef PyExperiment exp;
        for exp in self._exprmnts:
            vec.push_back((<PyExperiment>exp)._exprmnt)
        self._rdcntr = HammingReadCounter.getReadCounter(vec, mismatches, minlength, barcode_mismatches)

    @property
    def minimum_unique_barcode_length(self):
        return (<HammingReadCounter*>self._rdcntr).minimumUniqueBarcodeLength()

    @property
    def allowed_barcode_mismatches(self):
        return (<HammingReadCounter*>self._rdcntr).allowedBarcodeMismatches()

    @property
    def unique_forward_barcodes(self):
        return (<HammingReadCounter*>self._rdcntr).uniqueForwardBarcodes()

    @property
    def unique_reverse_barcodes(self):
        return (<HammingReadCounter*>self._rdcntr).uniqueReverseBarcodes()

cdef class PySeqlevReadCounter(PyReadCounter):
    def __cinit__(self, list experiments, int mismatches=1, int barcode_mismatches=0):
        cdef vector[Experiment*] vec
        cdef PyExperiment exp;
        for exp in self._exprmnts:
            vec.push_back((<PyExperiment>exp)._exprmnt)
        self._rdcntr = SeqlevReadCounter.getReadCounter(vec, mismatches, barcode_mismatches)

    @property
    def allowed_barcode_mismatches(self):
        return (<SeqlevReadCounter*>self._rdcntr).allowedBarcodeMismatches()

