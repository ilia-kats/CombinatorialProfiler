# cython: c_string_type=unicode, c_string_encoding=utf8

from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libc.stdint cimport *

import csv

cdef extern from "cReadCounter.h" nogil:
    cdef cppclass BarcodeSet:
        BarcodeSet(const unordered_map[string, vector[string]]&, const unordered_map[string, vector[string]]) except+

        const unordered_map[string, vector[string]] fw
        const unordered_map[string, vector[string]] rev

    cdef cppclass ReadCounter:
        ReadCounter(const string&, BarcodeSet*, BarcodeSet*) except+
        void countReads(const string&, const string&, int)
        const unordered_map[pair[string, string], unordered_map[string, uint64_t]]& getCounts()
        uint64_t read()
        uint64_t counted()
        uint64_t unmatchedInsert()
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
            codes = {}
            for r in reader:
                codes[r[0]] = r[1].upper().encode()
        pyfw = self._generate_unique(codes)
        pyrev = self._make_reverse(pyfw)

        cdef unordered_map[string, vector[string]] fw
        cdef unordered_map[string, vector[string]] rev
        for k,v in pyfw.items():
            fw[k.encode()] = v
        for k,v in pyrev.items():
            rev[k.encode()] = v
        self._bset = new BarcodeSet(fw, rev)

    def __dealloc__(self):
        del self._bset

    def _generate_unique(self, codes):
        fw = {}
        for k, b in codes.items():
            bcodes = []
            for s in range(0, len(b)):
                found = False
                for key, barcode in codes.items():
                    if not key == k and barcode.find(b[s:]) != -1:
                        found = True
                        break
                if not found:
                    bcodes.append(b[s:])
            fw[k] = tuple(bcodes)
        return fw

    def _make_reverse(self, fw):
        rev = {}
        for k, v in fw.items():
            rev[k] = tuple(PyBarcodeSet.reverse_compl(bcode) for bcode in v)
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
    cdef readonly unicode insertseq
    cdef readonly PyBarcodeSet fw
    cdef readonly PyBarcodeSet rev

    def __cinit__(self, unicode insertseq, PyBarcodeSet barcodes_fw, PyBarcodeSet barcodes_rev):
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
        self._rdcntr = new ReadCounter(insertseq.encode(), fw, rev)

    def __dealloc__(self):
        del self._rdcntr

    def countReads(self, unicode fpath, unicode unmatchedpattern, threads=1):
        self._rdcntr.countReads(fpath.encode(), unmatchedpattern.encode(), threads)

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

