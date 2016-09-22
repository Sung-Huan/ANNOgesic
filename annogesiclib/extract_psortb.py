import sys
from annogesiclib.gff3 import Gff3Parser


def import_psortb(seq_name, psortbs, local_name, local_score, type_, results):
    seq_datas = seq_name.split("__")
    seq_id = seq_datas[0]
    features = seq_datas[1].split("_")
    prot_id = "_".join(features[:-3])
    if type_ == "multi":
        psortbs.append({"seq_id": seq_id,
                        "protein_id": prot_id,
                        "strand": features[-3],
                        "start": int(features[-2]),
                        "end": int(features[-1]),
                        "local": "/".join(local_name),
                        "score": "/".join(local_score)})
    else:
        psortbs.append({"seq_id": seq_id,
                        "protein_id": prot_id,
                        "strand": features[-3],
                        "start": int(features[-2]),
                        "end": int(features[-1]),
                        "local": results[0],
                        "score": results[-1]})
    return {"datas": seq_datas, "features": features, "prot_id": prot_id}


def get_results(line, scores, psortbs, out_p, seq_name, fuzzy):
    '''print the results of psorb'''
    local_name = []
    local_score = []
    if len(line) == 0:
        pass
    elif "(This protein may have multiple localization sites.)" in line:
        results = line.split(" ")
        sort_scores = sorted(scores, key=lambda x: (x["score"]), reverse=True)
        first = True
        high_scores = []
        for score in sort_scores:
            if first:
                high_scores.append(score)
                first = False
                best_score = score
            else:
                if score["local"] != results[0]:
                    if score["score"] < (best_score["score"] - fuzzy):
                        break
                    else:
                        high_scores.append(score)
        for high_score in high_scores:
            local_name.append(high_score["local"])
            local_score.append(str(high_score["score"]))
        seq_datas = import_psortb(seq_name, psortbs, local_name,
                                  local_score, "multi", results)
        out_p.write("\t".join([seq_datas["datas"][0], seq_datas["prot_id"],
                               "\t".join(seq_datas["features"][-3:]),
                               "/".join(local_name),
                               "/".join(local_score)]) + "\n")
    else:
        results = line.split(" ")
        seq_datas = import_psortb(seq_name, psortbs, None,
                                  None, "unique", results)
        out_p.write("\t".join([seq_datas["datas"][0], seq_datas["prot_id"],
                               "\t".join(seq_datas["features"][-3:]),
                               results[0], results[-1]]) + "\n")
    return local_name, local_score


def get_information(psortb_table, out_p, fuzzy):
    '''get the information of psorb'''
    scores = []
    psortbs = []
    seq_name = None
    detects = {"score": False, "result": False}
    with open(psortb_table, "r") as p_h:
        for line in p_h:
            line = line.strip()
            if (line.startswith("--")) or \
               (line.startswith("Secondary localization(s):")):
                detects["result"] = False
            if detects["score"]:
                if "Final Prediction:" not in line:
                    datas = line.split(" ")
                    scores.append({"local": datas[0],
                                   "score": float(datas[-1])})
            if detects["result"]:
                local_name, local_score = get_results(
                    line, scores, psortbs, out_p, seq_name, fuzzy)
            if line.startswith("Final Prediction:"):
                detects["score"] = False
                detects["result"] = True
            if line.startswith("SeqID:"):
                seq_name = line.replace("SeqID: ", "")
                scores = []
            if line.startswith("Localization Scores:"):
                detects["score"] = True
    return psortbs


def print_gff(gffs, psortbs, out_m):
    for gff in gffs:
        detect = False
        for psortb in psortbs:
            if (gff.feature == "CDS") and \
               (gff.start == psortb["start"]) and \
               (gff.end == psortb["end"]) and \
               (gff.strand == psortb["strand"]):
                if "protein_id" in gff.attributes.keys():
                    if gff.attributes["protein_id"] == psortb["protein_id"]:
                        detect = True
                        break
                elif "locus_tag" in gff.attributes.keys():
                    if gff.attributes["locus_tag"] == psortb["protein_id"]:
                        detect = True
                        break
                else:
                    if gff.attributes["ID"] == psortb["protein_id"]:
                        detect = True
                        break
        if detect:
            gff.attribute_string = gff.attribute_string + \
                                   ";subcellular_localization=" + \
                                   psortb["local"]
            out_m.write("\t".join([gff.info_without_attributes,
                                   gff.attribute_string + "\n"]))
        else:
            out_m.write(gff.info + "\n")


def extract_psortb(psortb_table, out_psortb, merge_gff, out_merge, fuzzy):
    '''extract and re-generate the output information of psorb'''
    gffs = []
    if merge_gff:
        if out_merge is None:
            print("Error: Assign a name of output merged annotation file.")
            sys.exit()
        out_m = open(out_merge, "w")
        for entry in Gff3Parser().entries(open(merge_gff)):
            gffs.append(entry)
        gffs = sorted(gffs, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    else:
        out_m = None
    out_p = open(out_psortb, "w")
    out_p.write("#strain\tprotein\tstrand\tstart\tend\tlocation\tscore\n")
    psortbs = get_information(psortb_table, out_p, fuzzy)
    if merge_gff:
        print_gff(gffs, psortbs, out_m)
