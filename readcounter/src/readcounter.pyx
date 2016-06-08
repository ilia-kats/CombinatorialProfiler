# cython: c_string_type=unicode, c_string_encoding=utf8

from cython.operator cimport dereference as deref

from libcpp.unordered_map cimport unordered_map
from libcpp.unordered_set cimport unordered_set
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.memory cimport shared_ptr
from libc.stdint cimport *

import csv
from enum import Enum

import pandas as pd



cdef extern from "cReadCounter.h" nogil:
    ctypedef unordered_map[string, unordered_set[string]] InsertSet
    ctypedef unordered_map[string, unordered_map[string, vector[string]]] BarcodeSet

    ctypedef unordered_map[string, string] SequenceSet
    ctypedef unordered_map[string, unordered_map[string, double]] SortedCellCounts
    ctypedef unordered_map[pair[string, string], unordered_map[string, uint64_t]] Counts

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
        NDSIS ndsi
        SortedCellCounts sortedCells
        Counts counts

    cdef cppclass ReadCounter:
        ReadCounter(vector[Experiment*]) except+
        void countReads(const string&, const string&, int)
        uint64_t read()
        uint64_t unmatchedTotal()
        uint64_t unmatchedInsert()
        uint64_t unmatchedBarcodes()
        uint64_t unmatchedInsertSequence()
        uint64_t counted()
        uint64_t written()

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
            deref(self._exprmnt).insert = d['insert'].encode()
        else:
            raise RuntimeError("Insert sequence missing from experiment %s" % deref(self._exprmnt).name)
        if 'barcodes_fw' in d:
            for k, v in d['barcodes_fw'].items():
                deref(self._exprmnt).fwBarcodeSet[v.encode()] = k.encode()
        if 'barcodes_rev' in d:
            for k,v in d['barcodes_rev'].items():
                deref(self._exprmnt).revBarcodeSet[v.encode()] = k.encode()
        if 'ndsi' in d:
            haveNdsi = False
            if d['ndsi'] == 'forward':
                deref(self._exprmnt).ndsi = forward
                haveNdsi = True
            elif d['ndsi'] == 'reverse':
                deref(self._exprmnt).ndsi = reverse
                haveNdsi = True
            else:
                deref(self._exprmnt).ndsi = noNDSI
            if haveNdsi and 'sortedcells' in d:
                for fw, revs in d['sortedcells'].items():
                    for r, v in revs.items():
                        deref(self._exprmnt).sortedCells[fw.encode()][r.encode()] = v

    @property
    def name(self):
        return deref(self._exprmnt).name

    @property
    def insert(self):
        return deref(self._exprmnt).insert

    @property
    def forward_barcodes(self):
        return deref(self._exprmnt).fwBarcodeSet

    @property
    def reverse_barcodes(self):
        return deref(self._exprmnt).revBarcodeSet

    @property
    def ndsi(self):
        return NDSIS(deref(self._exprmnt).ndsi)

    @property
    def sorted_cells(self):
        return deref(self._exprmnt).sortedCells

    @property
    def counts(self):
        return deref(self._exprmnt).counts


cdef class PyReadCounter:
    cdef ReadCounter *_rdcntr
    cdef list _exprmnts

    def __cinit__(self, list experiments):
        cdef vector[Experiment*] vec
        cdef PyExperiment exp;
        self._exprmnts = experiments
        for exp in self._exprmnts:
            vec.push_back((<PyExperiment>exp)._exprmnt)
        self._rdcntr = new ReadCounter(vec)

    def __dealloc__(self):
        del self._rdcntr

    def countReads(self, unicode fpath, unicode unmatchedpattern, threads=1):
        self._rdcntr.countReads(fpath.encode(), unmatchedpattern.encode(), threads)

    #def asDataFrames(self):
        #frames = {}

        #countsdict = self.counts
        #for ins, bcodes in countsdict.items():
            #insert = []
            #barcode_fw = []
            #barcode_rev = []
            #sequence = []
            #counts = []
            #for codes, seqs in bcodes.items():
                #for seq, cnts in seqs.items():
                    #insert.append(ins)
                    #barcode_fw.append(codes[0])
                    #barcode_rev.append(codes[1])
                    #sequence.append(seq)
                    #counts.append(cnts)
            #df = pd.DataFrame()
            #df['insert'] = pd.Series(insert, dtype='category')
            #if self.fw is not None and ins in self.fw:
                #df['barcode_fw'] = pd.Categorical(barcode_fw, categories=self.fw[ins].keys())
            #else:
                #df['barcode_fw'] = ''
            #if self.rev is not None and ins in self.rev:
                #df['barcode_rev'] = pd.Categorical(barcode_rev, categories=self.rev[ins].keys())
            #else:
                #df['barcode_rev'] = ''
            #df['sequence'] = sequence

            #if self.namedInserts and ins in self.namedInserts:
                #df['named_insert'] = [self.namedInserts[ins][s] for s in sequence]

            #df['counts'] = counts
            #frames[ins] = df
        #return frames

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

