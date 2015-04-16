#!/usr/bin/python

import os        
import sys
import csv

class parser_splice(object):
    def parser(self, input_file):
        sh = open(input_file, "r");
        for row in csv.reader(sh, delimiter="\t"):
            yield assign_value(row)

class assign_value(object):
    def __init__(self, row):
        self.strain = row[0]
        self.start = int(row[1])
        self.end = int(row[2])
        self.splice = row[3]
        splice = row[3].split(":")
        self.supported_reads = int(splice[1])
        self.start_site_reads = int(splice[2])
        self.end_site_reads = int(splice[3])
        self.splice_type = splice[4]
        self.situation = splice[5]
        self.strand = row[5]
        self.info = ("\t".join(row))
    def __str__(self):
        return "%s %s %s %s %s" % (self.strain, self.start, self.end, self.splice, self.strand)
