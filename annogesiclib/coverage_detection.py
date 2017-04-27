import copy


def coverage_comparison(cover, cover_sets, poss, first, strand, cover_pos):
    '''Seaching the lowest and highest coverage'''
    if first:
        first = False
        cover_sets["high"] = cover
        cover_sets["low"] = cover
        poss["high"] = cover_pos
        poss["low"] = cover_pos
    else:
        if cover_sets["high"] < cover:
            cover_sets["high"] = cover
            poss["high"] = cover_pos
            poss["low"] = cover_pos
            cover_sets["low"] = cover
        if ((strand == "+") and (poss["low"] >= poss["high"])) or \
           ((strand == "-") and (poss["low"] <= poss["high"])):
            if cover_sets["low"] > cover:
                cover_sets["low"] = cover
                poss["low"] = cover_pos
        elif ((strand == "+") and (poss["low"] < poss["high"])) or \
             ((strand == "-") and (poss["low"] > poss["high"])):
            poss["low"] = cover_pos
            cover_sets["low"] = cover
    return first


def get_repmatch(replicates, cond):
    '''deal with the replicate match'''
    detect_all = False
    for rep in replicates:
        if "all" in rep:
            detect_all = True
            rep = int(rep.split("_")[-1])
            break
    if not detect_all:
        for match in replicates:
            if cond.split("_")[0] == match.split("_")[0]:
                rep = int(match.split("_")[-1])
    return rep


def define_cutoff(coverages, median, utr_type):
    '''get the cutoff'''
    cutoffs = {}
    if coverages[utr_type] == "mean":
        for track, values in median.items():
            cutoffs[track] = values["mean"]
    else:
        for track, values in median.items():
            cutoffs[track] = values["median"]
    return cutoffs


def check_notex(cover, texs, cutoff, notex):
    '''Check the cutoff of average coverage for TEX+ and TEX- libs'''
    if notex is not None:
        if len(texs) != 0:
            for keys in texs.keys():
                tracks = keys.split("@AND@")
                if cover["track"] == tracks[0]:
                    if cover["avg"] > cutoff:
                        return True
                elif cover["track"] == tracks[1]:
                    if cover["avg"] > notex:
                        return True
        else:
            if cover["avg"] > cutoff:
                return True
    else:
        if cover["avg"] > cutoff:
            return True


def run_tex(cover, texs, check_texs, tex_notex, type_,
            detect_num, poss, target_datas):
    '''Check the position of different libs'''
    if (cover["type"] == "tex") or (cover["type"] == "notex"):
        for key in texs.keys():
            if cover["track"] in key:
                texs[key] += 1
                check_texs[key].append(cover)
                if texs[key] >= tex_notex:
                    if type_ == "sRNA_utr_derived":
                        if detect_num == 0:
                            poss["start"] = cover["final_start"]
                            poss["end"] = cover["final_end"]
                        else:
                            exchange_start_end(poss, cover)
                    detect_num += 1
                    if cover not in target_datas:
                        target_datas.append(cover)
                    if tex_notex != 1:
                        if check_texs[key][0] not in target_datas:
                            target_datas.append(check_texs[key][0])
                            if type_ == "sRNA_utr_derived":
                                exchange_start_end(poss, check_texs[key][0])
    elif cover["type"] == "frag":
        if type_ == "sRNA_utr_derived":
            if detect_num == 0:
                poss["start"] = cover["final_start"]
                poss["end"] = cover["final_end"]
            else:
                exchange_start_end(poss, cover)
        detect_num += 1
        target_datas.append(cover)
    return detect_num


def check_tex(template_texs, covers, target_datas, notex, type_, poss, median,
              coverages, utr_type, cutoff_coverage, tex_notex):
    '''Check the candidates for TEX+/- libs 
    (should be detected in one or both of libs)'''
    detect_num = 0
    check_texs = {}
    texs = copy.deepcopy(template_texs)
    for key, num in texs.items():
        check_texs[key] = []
    for cover in covers:
        run_check_tex = False
        if type_ == "sRNA_utr_derived":
            cutoffs = define_cutoff(coverages, median, utr_type)
            if cover["track"] in cutoffs.keys():
                if cover["avg"] > cutoffs[cover["track"]]:
                    run_check_tex = True
            else:
                run_check_tex = True
        elif type_ == "sORF":
            if cover["avg"] > coverages[cover["track"]]:
                run_check_tex = True
        elif (type_ == "terminator"):
            run_check_tex = True
        elif (type_ == "normal"):
            run_check_tex = check_notex(cover, texs, cutoff_coverage,
                                        notex)
        else:
            if cover["avg"] > cutoff_coverage:
                run_check_tex = True
        if run_check_tex:
            detect_num = run_tex(cover, texs, check_texs, tex_notex,
                                 type_, detect_num, poss, target_datas)
    return detect_num


def exchange_start_end(poss, cover):
    '''modify the start and end point. get the long one'''
    if poss["start"] > cover["final_start"]:
        poss["start"] = cover["final_start"]
    if poss["end"] < cover["final_end"]:
        poss["end"] = cover["final_end"]


def replicate_comparison(args_srna, srna_covers, strand, type_, median,
                         coverages, utr_type, notex, cutoff_coverage, texs):
    '''Check the number of replicates which fit the cutoff in order to remove 
    the candidates which only can be detected in few replicates.'''
    srna_datas = {"best": 0, "high": 0, "low": 0, "start": -1,
                  "end": -1, "track": "", "detail": [], "conds": {}}
    tmp_poss = {"start": -1, "end": -1, "pos": -1,
                "all_start": [], "all_end": []}
    detect = False
    for cond, covers in srna_covers.items():
        detect_num = check_tex(
                texs, covers, srna_datas["detail"], notex, type_, tmp_poss,
                median, coverages, utr_type, cutoff_coverage,
                args_srna.tex_notex)
        if ("texnotex" in cond):
            tex_rep = get_repmatch(args_srna.replicates["tex"], cond)
            if detect_num >= tex_rep:
                detect = True
        elif ("frag" in cond):
            frag_rep = get_repmatch(args_srna.replicates["frag"], cond)
            if detect_num >= frag_rep:
                detect = True
        if detect:
            detect = False
            if type_ == "sRNA_utr_derived":
                tmp_poss["all_start"].append(tmp_poss["start"])
                tmp_poss["all_end"].append(tmp_poss["end"])
            else:
                if strand == "+":
                    sort_datas = sorted(srna_datas["detail"],
                                        key=lambda k: (k["pos"]))
                else:
                    sort_datas = sorted(srna_datas["detail"],
                                        key=lambda k: (k["pos"]), reverse=True)
                srna_datas["pos"] = sort_datas[-1]["pos"]
            sort_datas = sorted(srna_datas["detail"], key=lambda k: (k["avg"]))
            avg = sort_datas[-1]["avg"]
            srna_datas["conds"][cond] = str(detect_num)
            if (avg > srna_datas["best"]):
                srna_datas["high"] = sort_datas[-1]["high"]
                srna_datas["low"] = sort_datas[-1]["low"]
                srna_datas["best"] = avg
                srna_datas["track"] = sort_datas[-1]["track"]
    if type_ == "sRNA_utr_derived":
        if len(tmp_poss["all_start"]) != 0:
            srna_datas["start"] = min(tmp_poss["all_start"])
            srna_datas["end"] = max(tmp_poss["all_end"])
        else:
            srna_datas["start"] = -1
            srna_datas["end"] = -1
    return srna_datas
