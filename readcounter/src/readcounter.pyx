# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libc.stdint cimport *

import csv

import pandas as pd

cdef extern from "cReadCounter.h" nogil:
    cdef cppclass BarcodeSet:
        BarcodeSet(const unordered_map[string, unordered_map[string, vector[string]]]&, const unordered_map[string, unordered_map[string, vector[string]]]) except+

        const unordered_map[string, unordered_map[string, vector[string]]] fw
        const unordered_map[string, unordered_map[string, vector[string]]] rev

    cdef cppclass ReadCounter:
        ReadCounter(const unordered_map[string, string]&, BarcodeSet*, BarcodeSet*) except+
        void countReads(const string&, const string&, int)
        const unordered_map[string, unordered_map[pair[string, string], unordered_map[string, uint64_t]]]& getCounts()
        uint64_t read()
        uint64_t counted()
        uint64_t unmatchedInsert()
        uint64_t insertsWithoutBarcodes()
        uint64_t unmatchedBarcodeFw()
        uint64_t unmatchedBarcodeRev()


cdef class PyBarcodeSet:
    ttable = bytes.maketrans(b'ATGC', b'TACG')
    cdef BarcodeSet *_bset

    def __cinit__(self, fpath):
        with open(fpath) as bcodefile:
            sample = bcodefile.read(1024)
            bcodefile.seek(0)
            s = csv.Sniffer()
            dialect = s.sniff(sample)
            has_header = s.has_header(sample)
            reader = csv.reader(bcodefile, dialect)
            if has_header:
                next(reader)
            codes = dict()
            for r in reader:
                if len(r) > 2:
                    insname = r[2]
                else:
                    insname = ""
                if insname not in codes:
                    codes[insname] = {}
                codes[insname][r[0]] = r[1].upper().encode()
        pyfw = self._generate_unique(codes)
        pyrev = self._make_reverse(pyfw)

        cdef unordered_map[string, unordered_map[string, vector[string]]] fw
        cdef unordered_map[string, unordered_map[string, vector[string]]] rev

        for i, b in pyfw.items():
            for k,v in b.items():
                fw[i.encode()][k.encode()] = v
        for i,b in pyrev.items():
            for k,v in b.items():
                rev[i.encode()][k.encode()] = v
        self._bset = new BarcodeSet(fw, rev)

    def __dealloc__(self):
        del self._bset

    def _generate_unique(self, codes):
        fw = {}
        for i, c in codes.items():
            fw[i] = {}
            for k,b in c.items():
                bcodes = []
                for s in range(0, len(b)):
                    found = False
                    for key, barcode in codes[i].items():
                        if not key == k and barcode.find(b[s:]) != -1:
                            found = True
                            break
                    if not found:
                        bcodes.append(b[s:])
                fw[i][k] = tuple(bcodes)
        return fw

    def _make_reverse(self, fw):
        rev = {}
        for i, b in fw.items():
            rev[i] = {}
            for k,v in b.items():
                rev[i][k] = tuple(PyBarcodeSet.reverse_compl(bcode) for bcode in v)
        return rev

    @staticmethod
    def reverse_compl(sequence):
        return sequence.translate(PyBarcodeSet.ttable)[::-1]

    @property
    def fw(self):
        return self._bset.fw

    @property
    def rev(self):
        return self._bset.rev

cdef class PyReadCounter:
    cdef ReadCounter *_rdcntr
    cdef readonly dict insertseq
    cdef readonly PyBarcodeSet fw
    cdef readonly PyBarcodeSet rev

    def __cinit__(self, dict insertseq, PyBarcodeSet barcodes_fw, PyBarcodeSet barcodes_rev):
        self.insertseq = insertseq
        self.fw = barcodes_fw
        self.rev = barcodes_rev

        cdef BarcodeSet *fw
        cdef BarcodeSet *rev
        if barcodes_fw is None:
            fw = NULL
        else:
            fw = barcodes_fw._bset
        if barcodes_rev is None:
            rev = NULL
        else:
            rev = barcodes_rev._bset
        cdef unordered_map[string, string] seqs
        for k, v in insertseq.items():
            seqs[k.encode()] = v.encode()
        self._rdcntr = new ReadCounter(seqs, fw, rev)

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
            if len(countsdict) > 1:
                df['insert'] = insert
            df['barcode_fw'] = barcode_fw
            df['barcode_rev'] = barcode_rev
            df['sequence'] = sequence
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
    def inserts_without_barcodes(self):
        return self._rdcntr.insertsWithoutBarcodes()

    @property
    def unmatched_barcode_fw(self):
        return self._rdcntr.unmatchedBarcodeFw()

    @property
    def unmatched_barcode_rev(self):
        return self._rdcntr.unmatchedBarcodeRev()

