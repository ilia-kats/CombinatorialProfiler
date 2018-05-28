# -*- coding: utf-8 -*-
import csv
import string
from io import IOBase

import numpy as np
import pandas as pd

def get_csv_reader(f):
    if not isinstance(f, IOBase):
        f = open(f, newline=None)
    sample = f.read(1024)
    f.seek(0)
    s = csv.Sniffer()
    dialect = s.sniff(sample)
    has_header = s.has_header(sample)
    reader = csv.reader(f, dialect)
    return (reader, has_header, dialect)

def readBarcodes(f):
    reader, has_header, dialect = get_csv_reader(f)
    if has_header:
        next(reader)
    tr = str.maketrans('', '', 'ATCG')
    i = 0
    for r in reader:
        if len(r) < 2:
            raise RuntimeError("Not enough columns")
        if i == 0:
            if not r[1].upper().translate(tr):
                sc = 1
                nc = 0
            elif not r[0].upper().translate(tr):
                sc = 0
                nc = 1
            else:
                raise RuntimeError("No DNA sequence found in first two columns")
        i += 1
        yield (r[nc], r[sc].upper())

_unmatchedlabelsexception = RuntimeError("Barcode labels don't match sorted cells labels")

def readCellCounts(f, barcodes_fw=None, barcodes_rev=None, force_string=False):
    header, dialect = get_csv_reader(f)[1:3]
    if header:
        header = 0
    else:
        header = None
    df = pd.read_csv(f,
                     header=header,
                     sep=dialect.delimiter,
                     doublequote=dialect.doublequote,
                     escapechar=dialect.escapechar,
                     quotechar=dialect.quotechar,
                     skipinitialspace=dialect.skipinitialspace)
    df_toreturn = df
    if force_string:
        f.seek(0)
        df_toreturn = pd.read_csv(f,
                                  header=header,
                                  sep=dialect.delimiter,
                                  doublequote=dialect.doublequote,
                                  escapechar=dialect.escapechar,
                                  quotechar=dialect.quotechar,
                                  skipinitialspace=dialect.skipinitialspace,
                                  dtype=str)
    types = df.get_dtype_counts()
    if types['float64'] != df.columns.size - 1:
        valuecol = np.where(df.dtypes == 'float64')[0]
        if valuecol.size > 1:
            raise RuntimeError("Unrecognized sorted cells format: Multiple numeric columns")
        if valuecol[0] == 0:
            raise RuntimeError("Unrecognized sorted cells format: No barcode labels given")
        if barcodes_fw is not None and len(barcodes_fw) and barcodes_rev is not None and len(barcodes_rev):
            indexcols = [valuecol[0] - 2, valuecol[0] - 1]
        else:
            df[df.columns.size] = ''
            indexcols = [valuecol - 1, df.columns.size - 1]

        if barcodes_fw is not None and len(barcodes_fw) and not set(df.iloc[:,indexcols[0]].squeeze()).isdisjoint(barcodes_fw):
            if barcodes_rev is not None and len(barcodes_rev) and not set(df.iloc[:,indexcols[1]].squeeze()).isdisjoint(barcodes_rev):
                index = df.columns[indexcols].values.tolist()
            else:
                raise _unmatchedlabelsexception
        elif barcodes_rev is not None and len(barcodes_rev) and not set(df.iloc[:,indexcols[0]].squeeze()).isdisjoint(barcodes_rev):
            index = df.columns[list(reversed(indexcols))].values.tolist()
        else:
            raise _unmatchedlabelsexception
        toreturn = df_toreturn.set_index(index)
        toreturn.index.rename(['barcode_fw', 'barcode_rev'], inplace=True)
        toreturn = toreturn.loc[:,df.columns[valuecol]]
    else:
        df = df.set_index(df.columns[0])
        if barcodes_fw is not None and len(barcodes_fw):
            if not set(df.index).isdisjoint(barcodes_fw):
                toreturn = df_toreturn.stack()
            elif not set(df.columns).isdisjoint(barcodes_fw):
                toreturn = df_toreturn.T.stack()
            else:
                raise _unmatchedlabelsexception
        elif barcodes_rev is not None and len(barcodes_rev):
            if not set(df_toreturn.index).isdisjoint(barcodes_rev):
                toreturn = df_toreturn.T.stack()
            elif not set(df.columns).isdisjoint(barcodes_rev):
                toreturn = df_toreturn.stack()
            else:
                raise _unmatchedlabelsexception
        else:
            raise _unmatchedlabelsexception

        toreturn.index.rename(['barcode_fw', 'barcode_rev'], inplace=True)
    return toreturn

