from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp.list cimport list
from libcpp.utility cimport pair

import csv

cdef extern from "cReadCounter.h" nogil:
    cdef cppclass BarcodeSet:
        BarcodeSet(const unordered_map[string, list[string]]&, const unordered_map[string, list[string]]) except+

        const unordered_map[string, list[string]] fw
        const unordered_map[string, list[string]] rev


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

        cdef unordered_map[string, list[string]] fw
        cdef unordered_map[string, list[string]] rev
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

