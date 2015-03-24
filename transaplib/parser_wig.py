#!/usr/bin/python

import os        
import sys
import csv


class parser_wig(object):
    def parser(self, input_file, strand):
        track= ""
        strain = ""
        wigs = []
        with open(input_file, "r") as w_f:
            for line in w_f:
                line = line.strip()
                line = line.split(" ")
                if (line[0] == "variableStep"):
                    strain = line[1].split("=")
                    strain = strain[1]
                    pre_pos = 0
                    first = True
                if (line[0] == "track"):
                    track = line[2].split("=")
                    track = track[1].replace("\"", "")
                    pre_pos = 0
                    first = True
                if (line[0] != "track") and (line[0] != "variableStep"):
                    if int(line[0]) - 1 != pre_pos:
                        for pos in range(pre_pos + 1, int(line[0])):
                            yield assign_value(pos, 0, strand, strain, track)
                        pre_pos = int(line[0])
                        first = True
                    if (int(line[0]) - 1 == pre_pos) or (first):
                        pre_pos = int(line[0])
                        first = False
                        yield assign_value(line[0], line[1], strand, strain, track)

class assign_value(object):
    def __init__(self, pos, coverage, strand, strain, track):
        self.pos = int(pos)
        if strand == "+":
            self.coverage = float(coverage)
        else:
            self.coverage = -1 * float(coverage)
        self.strand = strand
        self.strain = strain
        self.track = track
            
    def __str__(self):
        return "%s %s %s %s %s" % (self.pos, self.coverage, self.strand, self.strain, self.track)
