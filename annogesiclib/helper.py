import os
import sys
import copy
import shutil
import csv
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from annogesiclib.gff3 import Gff3Parser


class Helper(object):
    '''For some small and regular modules for ANNOgesic'''

    def __init__(self):
        self.gff3parser = Gff3Parser()

    def feature_without_notgene(self, entry):
        if (entry.feature != "gene") and (
            entry.feature != "exon") and (
            entry.feature != "source") and (
            entry.feature != "region") and (
            entry.feature != "repeat_region") and (
            entry.feature != "transcript") and (
            entry.feature != "STS") and (
            entry.feature != "remark"):
            utr_markers = r'[^\'\ \-35]'
            for sub_f in entry.feature.lower().split("utr"):
                match = (re.search(utr_markers, sub_f))
                if match is None:
                    pass
                else:
                    return True
        return False

    def get_strand_name(self, strand):
        '''change the strand name to f/r'''
        name = ""
        if strand == "+":
            name = "f"
        else:
            name = "r"
        return name

    def _fix_break_line(self, tar, prefix):
        '''fix the break line at the end of file'''
        tmp_out = open("tmp_file", "w")
        first = True
        with open(tar) as fh:
            for line in fh:
                line = line.strip()
                if (prefix == ">"):
                    if (prefix in line) and (first):
                        first = False
                    elif (prefix in line) and (not first) and (
                          not line.startswith(prefix)):
                        line = line.replace(prefix, "\n" + prefix)
                else:
                    row = line.split("\t")
                    if (len(row) > 9):
                        for strain in prefix:
                            if strain in line:
                                line = line.replace(strain, "\n" + strain)
                                break
                tmp_out.write(line + "\n")
        tmp_out.close()
        os.remove(tar)
        shutil.move("tmp_file", tar)

    def merge_file(self, ref, tar):
        '''merge two files'''
        os.system(" ".join(["cat", ref, ">>", tar]))
        if tar.endswith(".fa"):
            self._fix_break_line(tar, ">")
        elif tar.endswith(".gff"):
            strains = []
            fh = open(ref, "r")
            for row in csv.reader(fh, delimiter='\t'):
                if row[0] not in strains:
                    strains.append(row[0])
            fh.close()
            self._fix_break_line(tar, strains)

    def merge_blast_out(self, ref, tar):
        tmp_out = tar + "_tmp"
        out = open(tmp_out, "w")
        file_num = 0
        for file_ in (ref, tar):
            start = False
            finish = False
            with open(file_) as fh:
                for line in fh:
                    check_line = line.strip()
                    if (check_line.startswith("Query")):
                        start = True
                    if (not start) and (not finish) and (file_num == 0):
                        out.write(line)
                    if start and (not check_line.startswith("Database")):
                        out.write(line)
                    if (check_line.startswith("Database")) and start:
                        if file_num == 1:
                            out.write(line)
                        else:
                            start = False
                            finish = True
            file_num += 1
        os.remove(tar)
        os.remove(ref)
        shutil.move(tmp_out, tar)

    def remove_all_content(self, folder, feature, data_type):
        '''remove all files in one folder'''
        for file_ in os.listdir(folder):
            remove = False
            if feature is None:
                remove = True
            else:
                if feature in file_:
                    remove = True
            if remove:
                target = os.path.join(folder, file_)
                if (data_type == "file") and os.path.isfile(target):
                    os.remove(target)
                elif (data_type == "dir") and os.path.isdir(target):
                    shutil.rmtree(target)

    def move_all_content(self, ref_folder, tar_folder, features):
        '''move all files form one folder to another one'''
        for file_ in os.listdir(ref_folder):
            move = False
            if (features is not None):
                move = True
                for feature in features:
                    if (feature not in file_):
                        move = False
            elif (features is None):
                move = True
            if move:
                shutil.move(os.path.join(ref_folder, file_),
                            os.path.join(tar_folder, file_))

    def remove_tmp(self, folder):
        '''remove the tmp folder'''
        if folder:
            if os.path.isdir(os.path.join(folder, "tmp")):
                shutil.rmtree(os.path.join(folder, "tmp"))
            self.remove_all_content(folder, "_folder", "dir")

    def remove_tmp_dir(self, folder):
        if folder is not None:
            if os.path.isdir(folder):
                shutil.rmtree(folder)

    def remove_tmp_dir_ori(self, ori_folder, folder, types):
        if folder is not None:
            if os.path.isdir(folder):
                for file_ in os.listdir(ori_folder):
                    for type_ in types:
                        ori_file = os.path.join(ori_folder, file_)
                        if (file_.endswith(type_)) and (
                                os.path.isfile(ori_file)):
                            shutil.move(ori_file, ori_file + "_old")
                for file_ in os.listdir(folder):
                    if os.path.isfile(os.path.join(folder, file_)):
                        shutil.move(os.path.join(folder, file_), ori_folder)
                shutil.rmtree(folder)

    def remove_wigs(self, wigs):
        '''remove the merge wig folder which is generated by ANNOgesic'''
        if wigs:
            folder = wigs.split("/")
            folder = "/".join(folder[:-1])
            if os.path.isdir(os.path.join(folder, "merge_wigs")):
                shutil.rmtree(os.path.join(folder, "merge_wigs"))
        self.remove_tmp(wigs)

    def get_correct_file(self, datas, feature, prefix, for_wig_type, libs):
        '''get the correct file by comparing the strain name'''
        detect = False
        for data in os.listdir(datas):
            if os.path.isfile(os.path.join(datas, data)):
                if for_wig_type is None:
                    if feature in data:
                        file_ = data[:-1 * len(feature)]
                        if prefix == file_:
                            detect = True
                            return os.path.join(datas, data)
                else:
                    filename = data.split("_STRAIN_")
                    if ("reverse" in data) and ("forward" in data):
                        print("Error: Assign reverse or forward wigs!!!")
                        sys.exit()
                    elif (prefix == filename[-1][:-1 * len(feature)]):
                        if (for_wig_type == "forward"):
                            for lib in libs:
                                if (filename[0] in lib) and (lib[-1] == "+"):
                                    return os.path.join(datas, data)
                        if (for_wig_type == "reverse"):
                            for lib in libs:
                                if (filename[0] in lib) and (lib[-1] == "-"):
                                    return os.path.join(datas, data)
        if detect:
            detect = False
        else:
            print("Warning: No proper file - " + prefix + feature)
            return None

    def check_make_folder(self, folder):
        '''make new folder (if the folder exists, 
        it will remove it and create new one)'''
        path = "/".join(folder.split("/")[:-1])
        gen_folder = folder.split("/")[-1]
        if gen_folder in os.listdir(path):
            shutil.rmtree(os.path.join(path, gen_folder))
        os.mkdir(os.path.join(path, gen_folder))

    def sort_gff(self, gff_file, out_file):
        gffs = []
        g_f = open(gff_file, "r")
        for entry in self.gff3parser.entries(g_f):
            gffs.append(entry)
        g_f.close()
        sort_gffs = sorted(gffs, key=lambda x: (x.seq_id, x.start,
                                                x.end, x.strand))
        out = open(out_file, "w")
        out.write("##gff-version 3\n")
        for gff in sort_gffs:
            out.write("\t".join([str(field) for field in [
                        gff.seq_id, gff.source, gff.feature, gff.start,
                        gff.end, gff.score, gff.strand, gff.phase,
                        gff.attribute_string]]) + "\n")
        out.close()

    def extract_gene(self, seq, start, end, strand):
        '''extract gene seqence'''
        fasta = ''
        if strand == "+":
            return seq[(int(start)-1):int(end)]
        else:
            rev_seq = seq[(int(start)-1):int(end)]
            fasta = self._reverse_seq(rev_seq)
            return fasta

    def _reverse_seq(self, rev_seq):
        '''deal with the reverse strand'''
        fasta = ""
        for base in rev_seq[::-1]:
            if base.upper() == 'A':
                fasta = fasta + 'T'
            elif base.upper() == 'T':
                fasta = fasta + 'A'
            elif base.upper() == 'C':
                fasta = fasta + 'G'
            elif base.upper() == 'G':
                fasta = fasta + 'C'
        return fasta

    def _add_element(self, list_, type_, gff):
        if type_ in gff.attributes.keys():
            list_.add(gff.attributes[type_])

    def check_uni_attributes(self, gff_file):
        '''check the attributes of gff filee, ie ID has to be unique'''
        print("Checking gff file of {0}".format(gff_file.split("/")[-1]))
        gffs = []
        fh = open(gff_file)
        lens = {}
        for entry in self.gff3parser.entries(fh):
            if (entry.feature == "source") or (
                    entry.feature == "region") or (
                    entry.feature == "remark"):
                if entry.seq_id not in lens.keys():
                    lens[entry.seq_id] = entry.end
                else:
                    if entry.end > lens[entry.seq_id]:
                        lens[entry.seq_id] = entry.end 
                    print("Warning: Detect repeated source/region/remark of {0}! "
                          "The longest end point is used as the length of the genome.".format(
                           entry.seq_id))
            gffs.append(entry)
        gffs = sorted(gffs, key=lambda x: (x.seq_id, x.start, x.end, x.strand))
        first = True
        ids = set()
        locus_tags = set()
        pre_gff = None
        for gff in gffs:
            if (gff.feature != "source") and (
                    gff.feature != "region") and (
                    gff.feature != "remark"):
                if gff.seq_id in lens.keys():
                    if (gff.end > lens[gff.seq_id]):
                        name = "".join([gff.feature, ":", str(gff.start), "-",
                                        str(gff.end), "_", gff.strand])
                        print("Error: The end point of " + name +
                              " is longer than the length of whole genome.")
                        print("Please check the gff files.")
                        sys.exit()
            if first:
                first = False
                self._add_element(ids, "ID", gff)
                self._add_element(locus_tags, "locus_tag", gff)
            else:
                if gff.seq_id == pre_gff.seq_id:
                    if "ID" in gff.attributes.keys():
                        if gff.attributes["ID"] in ids:
                            print("Warninng: Repeat ID {0} "
                                  "in gff file!!!".format(
                                      gff.attributes["ID"]))
                        else:
                            self._add_element(ids, "ID", gff)
                    if "locus_tag" in gff.attributes.keys():
                        if gff.attributes["locus_tag"] in ids:
                            print("Warning: Repeat locus_tag {0} "
                                  "in gff file!!!".format(
                                      gff.attributes["locus_tag"]))
                        else:
                            self._add_element(locus_tags, "locus_tag", gff)
            pre_gff = copy.copy(gff)
        fh.close()

    def _read_fasta(self, fasta_file):
        seq = ""
        with open(fasta_file, "r") as seq_f:
            for line in seq_f:
                line = line.strip()
                if line.startswith(">"):
                    continue
                else:
                    seq = seq + line
        return seq

    def get_seq(self, gff_file, fasta_file, out_file):
        gff_f = open(gff_file, "r")
        out = open(out_file, "w")
        seq = self._read_fasta(fasta_file)
        num = 0
        for entry in self.gff3parser.entries(gff_f):
            gene = self.extract_gene(seq, entry.start, entry.end, entry.strand)
            if "ID" in entry.attributes.keys():
                id_ = entry.attributes["ID"]
            else:
                id_ = entry.feature + str(num)
            out.write(">{0}|{1}|{2}|{3}|{4}\n{5}\n".format(
                      id_, entry.seq_id, entry.start,
                      entry.end, entry.strand, gene))
            num += 1
        gff_f.close()
        out.close()

    def get_cds_seq(self, gff_file, fasta_file, out_file):
        seq = self._read_fasta(fasta_file)
        out = open(out_file, "w")
        cdss = []
        gh = open(gff_file)
        for entry in self.gff3parser.entries(gh):
            if entry.feature == "CDS":
                cdss.append(entry)
        cdss = sorted(cdss, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
        for entry in cdss:
            cds = self.extract_gene(seq, entry.start, entry.end, entry.strand)
            if "protein_id" in entry.attributes.keys():
                protein_id = entry.attributes["protein_id"]
            elif "locus_tag" in entry.attributes.keys():
                protein_id = entry.attributes["locus_tag"]
            else:
                protein_id = entry.attributes["ID"]
            out.write("_".join([">" + entry.seq_id, "_" + protein_id,
                      entry.strand, str(entry.start), str(entry.end)]) + "\n")
            out.write(cds + "\n")
        out.close()
        gh.close()

    def translation(self, dna_file, protein_file):
        '''translate the DNA to residues'''
        out = open(protein_file, "w")
        with open(dna_file) as d_h:
            for seq in d_h:
                seq = seq.strip()
                if seq.startswith(">"):
                    out.write(seq + "\n")
                else:
                    dna = Seq(seq, generic_dna)
                    out.write(str(dna.translate()) + "\n")
        out.close()
