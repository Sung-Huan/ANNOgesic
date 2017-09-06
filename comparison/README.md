The scripts for comparison between ANNOgesic predictions and several databases
------------------------------------------------------------------------------

1. Please download the data from RegulonDB (http://regulondb.ccg.unam.mx/menu/download/datasets/index.jsp), 
EcoCyc (https://ecocyc.org/site-search.shtml) or DOOR2 (http://csbl.bmb.uga.edu/DOOR/index.php).

2. In order to make the comparison more reliable, please remove the non-expressed 
features of databases. Otherwise the performances will be influenced by the
non-expressed features.

3. For terminators, we also suggest to remove the terminators from the databases 
which does not contain a coverage significant decrease.

4. For sORF, we used the study of Hemm et. al (2010) as benchmarking set. Please 
convert the sORF information to Gff3 format.

5. The number of CRISPRs and riboswitches are few in databases, the manual 
comparison can be implemented easily. Thus, no scrips for the comparison is provided.
