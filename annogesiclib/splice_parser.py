import csv


class SpliceParser(object):
    '''parser the splice data of segemehl'''

    def parser(self, splice_fh):
        for row in csv.reader(splice_fh, delimiter="\t"):
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
        return "{0} {1} {2} {3} {4}".format(
               self.strain, self.start, self.end, self.splice, self.strand)
