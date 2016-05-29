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
        void countReads(const string&, int)
        const unordered_map[pair[string, string], unordered_map[string, uint64_t]]& getCounts()


cdef class PyBarcodeSet:
    ttable = bytes.maketrans(b'ATGC', b'TACG')
    cdef BarcodeSet *_bset
    cdef readonly dict fw
    cdef readonly dict rev

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
        self._generate_unique(codes)
        self._make_reverse()

        cdef unordered_map[string, vector[string]] fw
        cdef unordered_map[string, vector[string]] rev
        for k,v in self.fw.items():
            fw[k.encode()] = v
        for k,v in self.rev.items():
            rev[k.encode()] = v
        self._bset = new BarcodeSet(fw, rev)

    def __dealloc__(self):
        del self._bset

    def _generate_unique(self, codes):
        self.fw = {}
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
            self.fw[k] = tuple(bcodes)

    def _make_reverse(self):
        self.rev = {}
        for k, v in self.fw.items():
            self.rev[k] = tuple(PyBarcodeSet.reverse_compl(bcode) for bcode in v)

    @staticmethod
    def reverse_compl(sequence):
        return sequence.translate(PyBarcodeSet.ttable)[::-1]

cdef class PyReadCounter:
    cdef ReadCounter *_rdcntr
    cdef readonly string insertseq
    cdef readonly PyBarcodeSet fw
    cdef readonly PyBarcodeSet rev

    def __cinit__(self, string insertseq, PyBarcodeSet barcodes_fw, PyBarcodeSet barcodes_rev):
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
        self._rdcntr = new ReadCounter(insertseq, fw, rev)

    def __dealloc__(self):
        del self._rdcntr

    def countReads(self, fpath, threads=1):
        self._rdcntr.countReads(fpath, threads)

    @property
    def counts(self):
        return self._rdcntr.getCounts()

