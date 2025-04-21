#!/usr/bin/python

from intervaltree import Interval
from collections import defaultdict
from math import isnan
import logging

logger = logging.getLogger('diff_exon')

class Exons:
    """
    Class for exons and associated scores.
    Maps input from various file types to genome coordinates.
    """
    def __init__(self):
        self.dexon2scores = {}
        self.num_dexons = 0

    def read_dexons(self, exon_file, index, skip=0, dpsi_cut=0.05, sigscore_cut=0.05, include_nas=True):
        """
        exons of interest together with additional numerical info
        """
        logger.info("Reading in input exons...")
        with open(exon_file) as f:
            for i in range(skip):
                f.readline()
                
            for line in f:
                l = line.strip().split('\t')

                chrom = l[2]
                strand = l[3]
                # rmats uses 0-indexing
                if index == 0:
                    start = str(int(l[4]) + 1)
                else: # 1-indexing for suppa2
                    start = str(int(l[4]))
                end = l[5]

                dpsi = float('nan')
                try:
                    dpsi = float(l[8])
                except:
                    logger.debug("%s:%s-%s: dpsi %s not a number", chr, str(start), str(end), l[8])

                pval = float('nan')
                if len(l) > 9:
                    try:
                        pval = float(l[9])
                    except:
                        logger.debug("%s:%s-%s: pval %s not a number", chr, str(start), str(end), l[9])

                if abs(dpsi) >= dpsi_cut and (pval <= sigscore_cut or (isnan(pval) and include_nas)):
                    self.dexon2scores[(chrom, start, end, strand)] = (dpsi, pval)
        
        self.num_dexons = len(self.dexon2scores)
    
    def __len__(self):
        return self.num_dexons
