class WigParser(object):
    '''parser the wiggle file based on 
    strain, track, position and coverage'''

    def parser(self, wig_fh, strand):
        track = ""
        strain = ""
        for line in wig_fh.readlines():
            line = line.strip()
            if (len(line) != 0) and (not line.startswith("#")):
                datas = line.split(" ")
                if (datas[0] == "variableStep") or (datas[0] == "fixedStep"):
                    strain = datas[1].split("=")
                    strain = strain[1].strip()
                    pre_pos = 0
                    first = True
                if (datas[0] == "track"):
                    track = datas[2].split("=")
                    track = track[1].replace("\"", "")
                    pre_pos = 0
                    first = True
                if (datas[0] != "track") and (
                        datas[0] != "variableStep") and (
                        datas[0] != "fixedStep"):
                    if len(datas) != 2:
                        datas = line.split("\t")
                    if int(datas[0]) - 1 != pre_pos:
                        for pos in range(pre_pos + 1, int(datas[0])):
                            yield AssignValue(pos, 0, strand, strain, track)
                        pre_pos = int(datas[0])
                        first = True
                    if (int(datas[0]) - 1 == pre_pos) or (first):
                        pre_pos = int(datas[0])
                        first = False
                        yield AssignValue(datas[0], datas[1],
                                           strand, strain, track)


class AssignValue(object):

    def __init__(self, pos, coverage, strand, strain, track):
        self.pos = int(pos)
        if strand == "+":
            self.coverage = float(coverage)
        else:
            if float(coverage) < 0:
                self.coverage = -1 * float(coverage)
            else:
                self.coverage = float(coverage)
        self.strand = strand
        self.strain = strain
        self.track = track

    def __str__(self):
        return "{0} {1} {2} {3} {4}".format(
                self.pos, self.coverage, self.strand, self.strain, self.track)
