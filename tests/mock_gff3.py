class Create_generator(object):

    def __init__(self, gff, attributes, type_):
        if (type_ == "gff") or (type_ == "circ"):
            self.seq_id = gff["seq_id"]
            self.strain = gff["seq_id"]
            self.strand = gff["strand"]
            self.start = gff["start"]
            self.end = gff["end"]
            self.feature = gff["feature"]
            self.phase = gff["phase"]
            self.score = gff["score"]
            self.source = gff["source"]
            if type_ == "circ":
                self.supported_reads = gff["support"]
                self.start_site_reads = gff["start_site"]
                self.end_site_reads = gff["end_site"]
                self.situation = gff["situation"]
                self.splice_type = gff["splice_type"]
            self.attributes = {}
            for key, value in attributes.items():
                self.attributes[key] = value
            self.attribute_string = ";".join(
                ["=".join(items) for items in self.attributes.items()])
            self.info = "\t".join([str(field) for field in [
                            self.seq_id, self.source, self.feature, self.start,
                            self.end, self.score, self.strand, self.phase,
                            self.attribute_string]])
            self.info_without_attributes = "\t".join([str(field) for field in [
                            self.seq_id, self.source, self.feature, self.start,
                            self.end, self.score, self.strand, self.phase]])
        if type_ == "wig":
            self.coverage = gff["coverage"]
