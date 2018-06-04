from annogesiclib.gff3 import Gff3Parser


def get_proteins(datas, checks, blast_f, score_n):
    '''filter and import the protein hit of blast'''
    proteins = []
    nums = {"index": 0, "hypo": 0}
    for data in datas:
        if (nums["index"] % 4 == 0) and (nums["index"] != 0):
            if "[" in data:
                data1 = data.split("[")
                data2 = data1[1].split("]")
                strain = data2[0].strip()
            else:
                data1 = data.split("Length")
                strain = "NA"
            name = data1[0].strip()
            tag = datas[nums["index"] - 1]
            if ("hypothetical" in name) or (
                    "Hypothetical" in name) or (
                    "unknown" in name) or (
                    "Unknown" in name) or (
                    "predicted coding region" in name) or (
                    "Predicted coding region" in name) or (
                    "PREDICTED:" in name) or (
                    "putative" in name) or ("Putative" in name):
                nums["hypo"] += 1
            if not checks["detect"]:
                for line in blast_f:
                    line = line.strip()
                    if ("Expect" in line) and ("Score" in line) and (
                            "Method" in line):
                        e_value = line.split(",")[1].split(" ")[-1]
                        score = line.split("Score = ")[-1].split(" bits ")[0].strip()
                        if score_n is None:
                            checks["detect"] = True
                            checks["print"] = True
                            proteins.append({"name": name, "strain": strain,
                                             "e": e_value, "tag": tag,
                                             "score": score})
                        elif (score_n is not None) and (float(score) >= score_n):
                            checks["detect"] = True
                            checks["print"] = True
                            proteins.append({"name": name, "strain": strain,
                                             "e": e_value, "tag": tag,
                                             "score": score})
                        break
            else:
                if score_n is None:
                    proteins.append({"name": name, "strain": strain,
                                     "e": e_value, "tag": tag, "score": score})
                elif (score_n is not None) and (float(score) >= score_n):
                    proteins.append({"name": name, "strain": strain,
                                     "e": e_value, "tag": tag, "score": score})
        nums["index"] += 1
    return proteins, nums


def detect_hypo(proteins, blasts, type_):
    '''remove the hit which is hypothetical protein or unknown'''
    protein_names = {}
    for protein in proteins:
        name = protein["name"].replace("\n", "")
        if ("hypothetical" not in name) and (
                "Hypothetical" not in name) and (
                "Unknown" not in name) and (
                "unknown" not in name) and (
                "Predicted coding region" not in name) and (
                "predicted coding region" not in name) and (
                "PREDICTED:" not in name) and ("putative" not in name) and (
                "Putative" not in name):
            if (name not in protein_names.keys()):
                protein_names[name] = []
            protein_names[name].append(protein["tag"])
            if type_ != "equal":
                blasts["blast"] = True
    return protein_names, protein["e"], protein["score"]


def detect_nr(line, blast_f, out_t, blasts, prefix, score_n):
    '''detect the hit in nr database'''
    checks = {"print": False, "detect": False}
    if line.startswith(">"):
        info = line.replace(">", "")
        for line in blast_f:
            line = line.strip()
            if len(line) != 0:
                info = " ".join([info, line])
            else:
                break
        datas = info.split("|")
        proteins, nums = get_proteins(datas, checks, blast_f, score_n)
        if checks["print"]:
            if blasts["hit_num"] < 3:
                protein_names, e, score = detect_hypo(proteins, blasts, "low")
                if len(protein_names) != 0:
                    for key, value in protein_names.items():
                        out_t.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                                    prefix, key, ",".join(value), e, score))
                    blasts["hit_num"] += 1
    if checks["print"]:
        return 2
    else:
        return 1


def detect_srna(line, blast_f, out_t, blasts, prefix, score_s):
    '''detect hit in sRNA database'''
    print_ = False
    blasts["name"] = ""
    if line[0] == ">":
        name_complete = False
        blasts["name"] = line[1:]
        blasts["hit_num"] += 1
        for line in blast_f:
            line.strip()
            if line.startswith("Length="):
                name_complete = True
            if not name_complete:
                blasts["name"] = " ".join([blasts["name"], line])
            if "Expect =" in line:
                e_value = line.split(" ")[-1].strip()
                score = line.split("Score = ")[-1].split(" bits ")[0].strip()
                if score_s is None:
                    print_ = True
                elif (score_s is not None) and (float(score) >= score_s):
                    print_ = True
                break
    if print_:
        blasts["name"] = blasts["name"].lstrip().replace("\n", "")
        out_t.write("{0}\t{1}\t{2}\t{3}\n".format(
                    prefix, blasts["name"], e_value, score))
        blasts["blast"] = True
    if print_:
        return 2
    else:
        return 1


def read_gff(srna_file, data_type):
    srnas = []
    srna_f = open(srna_file, "r")
    for entry in Gff3Parser().entries(srna_f):
        attributes = {}
        for key, value in entry.attributes.items():
            if (data_type == "sRNA") and (
               "sRNA_hit" not in key):
                attributes[key] = value
            elif (data_type == "nr") and (
                  "nr_hit" not in key):
                attributes[key] = value
        entry.attributes = attributes
        attribute_string = ";".join(
            ["=".join(items) for items in entry.attributes.items()])
        entry.info = "\t".join([entry.info_without_attributes,
                                attribute_string])
        srnas.append(entry)
    srnas = sorted(srnas, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    return srnas


def print_file(database, out_f, info, srna_hit, nr_hit):
    if database == "sRNA":
        out_f.write("{0};sRNA_hit={1}\n".format(info, srna_hit))
    elif database == "nr":
        out_f.write("{0};nr_hit={1}\n".format(info, nr_hit))


def gen_out_flie(blasts, out_t, prefix, out_f, database, srna, names):
    if not blasts["blast"]:
        out_t.write("{0}\tNA\n".format(prefix))
        print_file(database, out_f,
                   srna.info, "NA", "NA")
    else:
        print_file(database, out_f, srna.info,
                   len(names), blasts["hit_num"])


def get_whole_query(line, blast_f):
    whole_line = line
    if (not line.endswith("+")) and (
            not line.endswith("-")):
        for line in blast_f:
            line = line.strip()
            whole_line = whole_line + line
            if (line.endswith("+")) or (
                    line.endswith("-")):
                return whole_line
    else:
        return whole_line


def extract_blast(blast_result, srna_file, output_file,
                  output_table, database, score_s, score_n):
    '''extract the result of blast'''
    out_f = open(output_file, "w")
    out_t = open(output_table, "w")
    out_f.write("##gff-version 3\n")
    srnas = read_gff(srna_file, database)
    for srna in srnas:
        blasts = {"hit_num": 0, "blast": False, "name": ""}
        names = []
        prefix = "\t".join([srna.seq_id, srna.attributes["ID"],
                           srna.strand, str(srna.start), str(srna.end)])
        print_ = 0
        with open(blast_result, "r") as blast_f:
            for line in blast_f:
                line = line.strip()
                if line.startswith("Query= "):
                    if print_ == 2:
                        print_ = 0
                    elif print_ == 1:
                        print_ = 0
                        out_t.write("{0}\tNA\n".format(prefix))
                    line = get_whole_query(line, blast_f)
                    go_out = False
                    query = line.split("=")[1].strip()
                    if (query == ("|".join([
                            srna.attributes["ID"], srna.seq_id,
                            str(srna.start), str(srna.end), srna.strand]))):
                        for line in blast_f:
                            line = line.strip()
                            if line.find("No hits found") != -1:
                                print_file(database, out_f,
                                           srna.info, "NA", "NA")
                                out_t.write("{0}\tNA\n".format(prefix))
                                break
                            elif line.find("Sequences producing "
                                           "significant alignments:") != -1:
                                for line in blast_f:
                                    line = line.strip()
                                    if len(line) != 0:
                                        if line.startswith(
                                                "Effective search space"):
                                            go_out = True
                                            break
                                        if database == "sRNA":
                                            p = detect_srna(
                                                line, blast_f, out_t,
                                                blasts, prefix, score_s)
                                            if p:
                                            	print_ = p
                                            if (len(blasts["name"]) > 0):
                                                if blasts["name"] not in names:
                                                    names.append(
                                                        blasts["name"])
                                        elif database == "nr":
                                            p = detect_nr(
                                                line, blast_f, out_t, blasts,
                                                prefix, score_n)
                                            if p:
                                                print_ = p
                                gen_out_flie(blasts, out_t, prefix, out_f,
                                             database, srna, names)
                                blasts["hit_num"] = 0
                                break
                        if go_out:
                            break
    out_f.close()
    out_t.close()


def extract_energy(srna_file, sec_file, out_file):
    '''extract the folding energy of sRNA'''
    s_f = open(srna_file, "r")
    check = False
    get_length = False
    out = open(out_file, "w")
    out.write("##gff-version 3\n")
    for srna in Gff3Parser().entries(s_f):
        with open(sec_file, "r") as d_f:
            for structure in d_f:
                structure = structure.rstrip('\n')
                if get_length:
                    length = len(structure)
                    get_length = False
                if (structure.startswith(">")):
                    if (("|".join([srna.attributes["ID"], srna.seq_id,
                         str(srna.start),
                         str(srna.end), srna.strand])) == structure[1:]) or (
                        ("|".join([srna.feature, srna.seq_id,
                         str(srna.start),
                         str(srna.end), srna.strand])) == structure[1:]):
                        check = True
                        get_length = True
                if (check) and (
                    (structure[0] == "(") or (
                     structure[0] == ")") or (
                     structure[0] == ".")) and (
                     structure[-1] == ")"):
                    check = False
                    data = structure.split(" ")
                    if (data[-1].find("(") == -1):
                        energy = float(data[-1][0:-1])
                    else:
                        energy = float(data[-1][1:-1])
                    out.write("{0};2d_energy={1:.4f}\n".format(
                              srna.info, (energy / float(length))))
                    break
    s_f.close()
    out.close()
