#!/usr/bin/env python
# module docstring

# Import section
import os, sys, pickle
from urllib.parse import unquote
from collections import Counter, namedtuple

# Globals/constants section

# Objects
GffEntry = namedtuple('GffEntry', (
    'seqid',
    'source',
    'gff_type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'attributes'
    )
)


class Genome:
    
    def __init__(self, header):
        self.chromosomes = {}
        self._parse_header(header)
        
    @property
    def n_chromosomes(self):
        return len(self.chromosomes)
    
    @property
    def n_exons(self):
        exon_count = 0
        for chrom in self.chromosomes.values():
            exon_count += chrom.n_exons
        return exon_count
    
    @property
    def length(self):
        length = 0
        for chrom in self.chromosomes.values():
            length += chrom.length
        return length
    
    def _parse_header(self, header):
        for metadata in header:
            tag, val = metadata[2:].split(maxsplit=1)
            tag = tag.replace('-', '_')
            setattr(self, tag, val)


class BaseGFF:
    
    def __init__(
        self,
        seqid,
        source,
        gfftype,
        start,
        end,
        score,
        strand,
        phase,
        attributes
                ):
        self.seqid = seqid
        self.source = source
        self.gfftype = gfftype
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
        self._init_feature()
    
    @property
    def length(self):
        return self.end - self.start + 1
    
    def _init_feature(self):
        pass
    
    def __str__(self):
        return '\t'.join([
            self.seqid,
            self.source,
            self.gfftype,
            str(self.start),
            str(self.end),
            self.score,
            self.strand,
            self.phase,
            str(self.attributes)
        ])
    
    def __repr__(self):
        return f'{type(self).__name__}({self.seqid}, {self.source}, {self.gfftype}, {self.start}, {self.end}, {self.score}, {self.strand}, {self.phase}, {self.attributes})'
# Functions

# operational code

# runtime behavior
if __name__ == '__main__':
    args = sys.argv
    infile = args[1]
    outfile_name = args[2]
    
    try:
        first_n = args[3]
    except IndexError:
        first_n = None

    genome = process_gff(infile, first_n = first_n)
    #TODO: pickle the genome
    