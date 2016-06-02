# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref

from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libc.stdint cimport *

import csv

import pandas as pd

cdef extern from "cReadCounter.h" nogil:
    ctypedef unordered_map[string, unordered_set[string]] InsertSet
    ctypedef unordered_map[string, unordered_map[string, vector[string]]] BarcodeSet

    cdef cppclass ReadCounter:
        ReadCounter(const unordered_map[string, string]&, BarcodeSet*, BarcodeSet*, InsertSet*) except+
        void countReads(const string&, const string&, int)
        const unordered_map[string, unordered_map[pair[string, string], unordered_map[string, uint64_t]]]& getCounts()
        uint64_t read()
        uint64_t counted()
        uint64_t unmatchedInsert()
        uint64_t unmatchedBarcodeFw()
        uint64_t unmatchedBarcodeRev()
        uint64_t unmatchedTotal()
        uint64_t written()

cdef dictToBarcodeSet(dict bset, BarcodeSet *out):
    cdef string istring
    cdef string kstring
    for i, c in bset.items():
        istring = i.encode()
        for k, v in c.items():
            kstring = k.encode()
            for b in v:
                deref(out)[istring][kstring].push_back(b.encode())

cdef dictToInsertSet(dict namedInserts, InsertSet *out):
    for i, v in namedInserts.items():
        for s, n in v.items():
            deref(out)[i.encode()].insert(s.encode())

cdef class PyReadCounter:
    cdef ReadCounter *_rdcntr
    cdef BarcodeSet _fw
    cdef BarcodeSet _rev
    cdef InsertSet _ins

    cdef readonly dict insertseq
    cdef readonly dict fw
    cdef readonly dict rev
    cdef readonly dict namedInserts

    def __cinit__(self, dict insertseq, dict barcodes_fw, dict barcodes_rev, dict insertsToMap):
        self.insertseq = insertseq
        self.fw = barcodes_fw
        self.rev = barcodes_rev
        self.namedInserts = insertsToMap

        cdef BarcodeSet *fw
        cdef BarcodeSet *rev
        cdef InsertSet *insset
        if barcodes_fw is None:
            fw = NULL
        else:
            fw = &self._fw
            dictToBarcodeSet(barcodes_fw, fw)
        if barcodes_rev is None:
            rev = NULL
        else:
            rev = &self._rev
            dictToBarcodeSet(barcodes_rev, rev)
        if insertsToMap is None:
            insset = NULL
        else:
            insset = &self._ins
            dictToInsertSet(insertsToMap, insset)

        cdef unordered_map[string, string] seqs
        for k, v in insertseq.items():
            seqs[k.encode()] = v.encode()
        self._rdcntr = new ReadCounter(seqs, fw, rev, insset)

    def __dealloc__(self):
        del self._rdcntr

    def countReads(self, unicode fpath, unicode unmatchedpattern, threads=1):
        self._rdcntr.countReads(fpath.encode(), unmatchedpattern.encode(), threads)

    def asDataFrames(self):
        frames = {}

        countsdict = self.counts
        for ins, bcodes in countsdict.items():
            insert = []
            barcode_fw = []
            barcode_rev = []
            sequence = []
            counts = []
            for codes, seqs in bcodes.items():
                for seq, cnts in seqs.items():
                    insert.append(ins)
                    barcode_fw.append(codes[0])
                    barcode_rev.append(codes[1])
                    sequence.append(seq)
                    counts.append(cnts)
            df = pd.DataFrame()
            df['insert'] = pd.Series(insert, dtype='category')
            df['barcode_fw'] = pd.Series(barcode_fw, dtype='category')
            df['barcode_rev'] = pd.Series(barcode_rev, dtype='category')
            df['sequence'] = sequence

            if self.namedInserts and ins in self.namedInserts:
                df['named_insert'] = [self.namedInserts[ins][s] for s in sequence]

            df['counts'] = counts
            frames[ins] = df
        return frames

    @property
    def counts(self):
        return self._rdcntr.getCounts()

    @property
    def read(self):
        return self._rdcntr.read()

    @property
    def counted(self):
        return self._rdcntr.counted()

    @property
    def unmatched_insert(self):
        return self._rdcntr.unmatchedInsert()

    @property
    def unmatched_barcode_fw(self):
        return self._rdcntr.unmatchedBarcodeFw()

    @property
    def unmatched_barcode_rev(self):
        return self._rdcntr.unmatchedBarcodeRev()

    @property
    def unmatched_total(self):
        return self._rdcntr.unmatchedTotal()

    @property
    def written(self):
        return self._rdcntr.written()

