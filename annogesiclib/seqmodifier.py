class SeqModifier(object):
    """Help to apply SNPs, insertion and deletions to a sequence."""

    def __init__(self, seq):
        self._seq = seq
        self._init_pos_dict()

    def seq(self):
        return self._seq

    def _init_pos_dict(self):
        self._org_pos_to_internal_pos = dict(
            [(pos, pos-1)
             for pos in range(1, len(self._seq) + 1)])

    def replace(self, pos, nucleotide):
        seq_as_list = list(self._seq)
        seq_as_list[self._org_pos_to_internal_pos[pos]] = nucleotide
        self._seq = "".join(seq_as_list)

    def remove(self, pos, num):
        int_pos = self._org_pos_to_internal_pos[pos]
        self._seq = self._seq[:int_pos] + self._seq[int_pos+1:]
        del(self._org_pos_to_internal_pos[pos])
        for pos in range(pos, len(self._seq) + 2):
            try:
                self._org_pos_to_internal_pos[pos] = (
                    self._org_pos_to_internal_pos[pos] - num)
            except KeyError:
                pass

    def insert(self, pos, nucleotide):
        """Insert after nucleotide of the given position"""
        int_pos = self._org_pos_to_internal_pos[pos]
        self._seq = self._seq[:int_pos] + nucleotide + self._seq[int_pos:]
        for pos in range(pos + 1, len(self._seq) + 1):
            try:
                self._org_pos_to_internal_pos[pos] = (
                    self._org_pos_to_internal_pos[pos] + len(nucleotide))
            except KeyError:
                pass

    def get_nucl(self, pos):
        return self._seq[self._org_pos_to_internal_pos[pos]]
