import sys
import os
import unittest
import shutil
from io import StringIO
sys.path.append(".")
import annogesiclib.get_Rfam_ribo as grr


class TestGetRfamRibo(unittest.TestCase):

    def setUp(self):
        self.example = Example()
        self.test_folder = "test_folder"
        if (not os.path.exists(self.test_folder)):
            os.mkdir(self.test_folder)

    def tearDown(self):
        if os.path.exists(self.test_folder):
            shutil.rmtree(self.test_folder)

    def test_rbs_from_rfam(self):
        out_file = os.path.join(self.test_folder, "test.out")
        ribo_table = os.path.join(self.test_folder, "ribo_table")
        with open(ribo_table, "w") as fh:
            fh.write(self.example.ribo_table)
        rfam_file = os.path.join(self.test_folder, "rfam")
        with open(rfam_file, "w") as fh:
            fh.write(self.example.rfam_file)
        grr.rbs_from_rfam(ribo_table, rfam_file, out_file)
        datas = []
        with open(out_file) as fh:
            for line in fh:
                line = line.strip()
                datas.append(line)
        self.assertEqual(set(datas), set(self.example.out.split("\n")))

class Example(object):
    ribo_table = """RF00162	SAM	SAM riboswitch box leader
RF00174	Cobalamin	Cobalamin riboswitch
RF00634	SAM-IV	S adenosyl methionine SAM riboswitch"""
    rfam_file = """INFERNAL1/a [1.1 | October 2013]
NAME     SAM
ACC      RF00162
STATES   338
NODES    91
CLEN     108
W        272
ALPH     RNA
RF       no
CONS     yes
MAP      yes
DATE     Thu Feb 20 18:18:58 2014
COM      [1] cmbuild -F CM SEED
COM      [2] cmcalibrate --mpi CM
PBEGIN   0.05
PEND     0.05
WBETA    1e-07
QDBBETA1 1e-07
QDBBETA2 1e-15
N2OMEGA  1.52588e-05
N3OMEGA  1.52588e-05
ELSELF   -0.08926734
//

INFERNAL1/a [1.1rc4 | June 2013]
NAME     SAM-IV
ACC      RF00634
STATES   362
NODES    95
CLEN     116
W        149
ALPH     RNA
RF       no
CONS     yes
MAP      yes
DATE     Tue Aug 20 22:41:06 2013
COM      [1] cmbuild -F CM SEED
COM      [2] cmcalibrate --mpi CM
PBEGIN   0.05
PEND     0.05
WBETA    1e-07
QDBBETA1 1e-07
QDBBETA2 1e-15
N2OMEGA  1.52588e-05
N3OMEGA  1.52588e-05
ELSELF   -0.08926734
NSEQ     40
EFFN     1.655273
CKSUM    531515698
NULL     0.000  0.000  0.000  0.000
GA       38.00
TC       38.60
NC       37.10
EFP7GF   -7.0145 0.71835
//

INFERNAL1/a [1.1rc4 | June 2013]
NAME     test
ACC      RF001111
STATES   362
NODES    95
CLEN     116
W        149
ALPH     RNA
RF       no
CONS     yes
MAP      yes
DATE     Tue Aug 20 22:41:06 2013
COM      [1] cmbuild -F CM SEED
COM      [2] cmcalibrate --mpi CM
PBEGIN   0.05
PEND     0.05
WBETA    1e-07
QDBBETA1 1e-07
QDBBETA2 1e-15
N2OMEGA  1.52588e-05
N3OMEGA  1.52588e-05
ELSELF   -0.08926734
NSEQ     40
EFFN     1.655273
CKSUM    531515698
NULL     0.000  0.000  0.000  0.000
GA       38.00
TC       38.60
NC       37.10
EFP7GF   -7.0145 0.71835

"""
    out = """INFERNAL1/a [1.1 | October 2013]
NAME     SAM
ACC      RF00162
STATES   338
NODES    91
CLEN     108
W        272
ALPH     RNA
RF       no
CONS     yes
MAP      yes
DATE     Thu Feb 20 18:18:58 2014
COM      [1] cmbuild -F CM SEED
COM      [2] cmcalibrate --mpi CM
PBEGIN   0.05
PEND     0.05
WBETA    1e-07
QDBBETA1 1e-07
QDBBETA2 1e-15
N2OMEGA  1.52588e-05
N3OMEGA  1.52588e-05
ELSELF   -0.08926734
//

INFERNAL1/a [1.1rc4 | June 2013]
NAME     SAM-IV
ACC      RF00634
STATES   362
NODES    95
CLEN     116
W        149
ALPH     RNA
RF       no
CONS     yes
MAP      yes
DATE     Tue Aug 20 22:41:06 2013
COM      [1] cmbuild -F CM SEED
COM      [2] cmcalibrate --mpi CM
PBEGIN   0.05
PEND     0.05
WBETA    1e-07
QDBBETA1 1e-07
QDBBETA2 1e-15
N2OMEGA  1.52588e-05
N3OMEGA  1.52588e-05
ELSELF   -0.08926734
NSEQ     40
EFFN     1.655273
CKSUM    531515698
NULL     0.000  0.000  0.000  0.000
GA       38.00
TC       38.60
NC       37.10
EFP7GF   -7.0145 0.71835
//"""

if __name__ == "__main__":
    unittest.main()

