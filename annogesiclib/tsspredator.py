import os
import sys
import shutil
from subprocess import call
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.converter import Converter
from annogesiclib.merge_manual import merge_manual_predict_tss
from annogesiclib.stat_TSSpredator import stat_tsspredator
from annogesiclib.plot_TSS_venn import plot_venn
from annogesiclib.validate_gene import validate_gff
from annogesiclib.stat_TA_comparison import stat_ta_tss
from annogesiclib.check_orphan import check_orphan
from annogesiclib.filter_TSS_pro import filter_tss_pro
from annogesiclib.filter_low_expression import filter_low_expression


class TSSpredator(object):

    def __init__(self, args_tss):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.master = os.path.join(args_tss.out_folder, "MasterTables")
        self.tmps = {"tss": "tmp_TSS", "ta_tss": "tmp_ta_tss", "tss_ta":
                     "tmp_tss", "tmp": "tmp"}
        if args_tss.ta_files is not None:
            self.tmps["ta"] = os.path.join(args_tss.ta_files, "tmp")
        else:
            self.tmps["ta"] = None
        self.gff_path = os.path.join(args_tss.gffs, "tmp")
        self.wig_path = os.path.join(args_tss.wig_folder, "tmp")
        self.fasta_path = os.path.join(args_tss.fastas, "tmp")
        self.stat_outfolder = os.path.join(args_tss.out_folder, "statistics")
        self.gff_outfolder = os.path.join(args_tss.out_folder, "gffs")

    def _assign_dict(self, lib_datas):
        return {"wig": lib_datas[0],
                "tex": lib_datas[1],
                "condition": int(lib_datas[2]),
                "replicate": lib_datas[3],
                "strand": lib_datas[4]}

    def _print_lib(self, lib_num, lib_list, out, wig_folder, prefix, rep_set):
        for num_id in range(1, lib_num+1):
            cond_list = []
            for lib in lib_list:
                if num_id == lib["condition"]:
                    cond_list.append(lib)
            cond_sort_list = sorted(cond_list, key=lambda k: k['replicate'])
            reps = []
            for cond in cond_sort_list:
                out.write("{0}_{1}{2} = {3}\n".format(
                          prefix, cond["condition"], cond["replicate"],
                          os.path.join(wig_folder, cond["wig"])))
                reps.append(cond["replicate"])
            for rep in sorted(rep_set):
                if rep not in reps:
                    out.write("{0}_{1}{2} = \n".format(
                              prefix, cond["condition"], rep))

    def _start_to_run(self, tsspredator_path, config_file, out_path, prefix):
        print("Running TSSpredator for " + prefix)
        out = open(os.path.join(out_path, "log.txt"), "w")
        err = open(os.path.join(out_path, "err.txt"), "w")
        call(["java", "-jar", tsspredator_path,
              config_file], stdout=out, stderr=err)
        out.close()
        err.close()

    def _import_lib(self, libs, wig_folder, project_strain_name,
                    out, gff, program, fasta):
        lib_dict = {"fp": [], "fm": [], "nm": [], "np": []}
        lib_num = 0
        rep_set = set()
        list_num_id = []
        print("Runniun {0} now...".format(program))
        for lib in libs:
            lib_datas = lib.split(":")
            if not lib_datas[0].endswith(".wig"):
                print("Error:Exist a not proper wig files!!")
                sys.exit()
            for wig in os.listdir(wig_folder):
                filename = wig.split("_STRAIN_")
                if (filename[0] == lib_datas[0][:-4]) and (
                        filename[1][:-4] == project_strain_name):
                    lib_datas[0] = wig
            if int(lib_datas[2]) > lib_num:
                lib_num = int(lib_datas[2])
            if lib_datas[3] not in rep_set:
                rep_set.add(lib_datas[3])
            if (lib_datas[1] == "tex") and (lib_datas[4] == "+"):
                lib_dict["fp"].append(self._assign_dict(lib_datas))
            elif (lib_datas[1] == "tex") and (lib_datas[4] == "-"):
                lib_dict["fm"].append(self._assign_dict(lib_datas))
            elif (lib_datas[1] == "notex") and (lib_datas[4] == "+"):
                lib_dict["np"].append(self._assign_dict(lib_datas))
            elif (lib_datas[1] == "notex") and (lib_datas[4] == "-"):
                lib_dict["nm"].append(self._assign_dict(lib_datas))
        for num_id in range(1, lib_num+1):
            out.write("annotation_{0} = {1}\n".format(num_id, gff))
        if program.lower() == "tss":
            self._print_lib(lib_num, lib_dict["fm"], out,
                            wig_folder, "fivePrimeMinus", rep_set)
            self._print_lib(lib_num, lib_dict["fp"], out,
                            wig_folder, "fivePrimePlus", rep_set)
        elif program.lower() == "processing_site":
            self._print_lib(lib_num, lib_dict["nm"], out,
                            wig_folder, "fivePrimeMinus", rep_set)
            self._print_lib(lib_num, lib_dict["np"], out,
                            wig_folder, "fivePrimePlus", rep_set)
        else:
            print("Error: Wrong program name!!!")
            sys.exit()
        for num_id in range(1, lib_num+1):
            out.write("genome_{0} = {1}\n".format(num_id, fasta))
        for num_id in range(1, lib_num+1):
            list_num_id.append(str(num_id))
        return lib_num, num_id, rep_set, lib_dict, list_num_id

    def _print_repmatch(self, args_tss, out):
        '''check replicate match'''
        if "all" in args_tss.repmatch:
            match = args_tss.repmatch.split("_")[-1]
            out.write("minNumRepMatches = {0}\n".format(match))
        else:
            nums = {}
            matchs = {}
            for match in args_tss.repmatch.split(","):
                lib = match.split("_")[0]
                rep = match.split("_")[-1]
                matchs[lib] = rep
                if rep not in nums.keys():
                    nums[rep] = 1
                else:
                    nums[rep] += 1
            for rep, num in nums.items():
                if num == max(nums.values()):
                    out.write("minNumRepMatches = {0}\n".format(rep))
                    max_rep = rep
                    break
            for lib, rep in matchs.items():
                if rep != max_rep:
                    out.write("minNumRepMatches_{0} = {1}\n".format(
                        lib, rep))

    def _gen_config(self, project_strain_name, args_tss, gff,
                    wig_folder, fasta, config_file):
        '''generation of config files'''
        master_folder = "MasterTable_" + project_strain_name
        out_path = os.path.join(self.master, master_folder)
        self.helper.check_make_folder(out_path)
        out = open(config_file, "w")
        out.write("TSSinClusterSelectionMethod = HIGHEST\n")
        out.write("allowedCompareShift = 1\n")
        out.write("allowedRepCompareShift = 1\n")
        lib_num, num_id, rep_set, lib_dict, list_num_id = \
            self._import_lib(args_tss.libs, wig_folder, project_strain_name,
                             out, gff, args_tss.program, fasta)
        out.write("idList = ")
        out.write(",".join(list_num_id) + "\n")
        out.write("maxASutrLength = 100\n")
        out.write("maxGapLengthInGene = 500\n")
        out.write("maxNormalTo5primeFactor = {0}\n".format(
                  args_tss.processing_factor))
        out.write("maxTSSinClusterDistance = {0}\n".format(
                  args_tss.cluster + 1))
        out.write("maxUTRlength = {0}\n".format(args_tss.utr_length))
        out.write("min5primeToNormalFactor = {0}\n".format(
                  args_tss.enrichment_factor))
        out.write("minCliffFactor = {0}\n".format(args_tss.factor))
        out.write("minCliffFactorDiscount = {0}\n".format(
                  args_tss.factor_reduction))
        out.write("minCliffHeight = {0}\n".format(args_tss.height))
        out.write("minCliffHeightDiscount = {0}\n".format(
                  args_tss.height_reduction))
        out.write("minNormalHeight = {0}\n".format(args_tss.base_height))
        self._print_repmatch(args_tss, out)
        out.write("minPlateauLength = 0\n")
        out.write("mode = cond\n")
        out.write("normPercentile = 0.9\n")
        if args_tss.program.lower() == "tss":
            self._print_lib(lib_num, lib_dict["nm"], out,
                            wig_folder, "normalMinus", rep_set)
            self._print_lib(lib_num, lib_dict["np"], out,
                            wig_folder, "normalPlus", rep_set)
        else:
            self._print_lib(lib_num, lib_dict["fm"], out,
                            wig_folder, "normalMinus", rep_set)
            self._print_lib(lib_num, lib_dict["fp"], out,
                            wig_folder, "normalPlus", rep_set)
        out.write("numReplicates = {0}\n".format(len(rep_set)))
        out.write("numberOfDatasets = {0}\n".format(lib_num))
        out.write("outputDirectory = {0}\n".format(out_path))
        for prefix_id in range(len(args_tss.output_prefixs)):
            out.write("outputPrefix_{0} = {1}\n".format(
                      prefix_id + 1, args_tss.output_prefixs[prefix_id]))
        out.write("projectName = {0}\n".format(project_strain_name))
        out.write("superGraphCompatibility = igb\n")
        out.write("texNormPercentile = 0.5\n")
        out.write("writeGraphs = 0\n")
        out.write("writeNocornacFiles = 0\n")
        out.close()

    def _convert_gff(self, prefixs, args_tss):
        for prefix in prefixs:
            out_file = os.path.join(self.gff_outfolder, "_".join([
                           prefix, args_tss.program]) + ".gff")
            gff_f = open(out_file, "w")
            out_path = os.path.join(self.master, "_".join([
                           "MasterTable", prefix]))
            if "MasterTable.tsv" not in os.listdir(out_path):
                print("Error:there is not MasterTable file in {0}".format(
                      out_path))
                print("Please check configuration file.")
            else:
                if args_tss.program.lower() == "processing":
                    feature = "processing_site"
                elif args_tss.program.lower() == "tss":
                    feature = "TSS"
                self.converter.convert_mastertable2gff(
                    os.path.join(out_path, "MasterTable.tsv"),
                    "ANNOgesic", feature, prefix, out_file)
            gff_f.close()

    def _merge_manual(self, tsss, args_tss):
        '''if manual detected TSS is provided, it can merge manual detected TSS 
        and TSSpredator predicted TSS'''
        self.helper.check_make_folder(os.path.join(os.getcwd(),
                                      self.tmps["tss"]))
        for tss in tsss:
            for gff in os.listdir(args_tss.gffs):
                if (gff[:-4] == tss) and (".gff" in gff):
                    break
            filename = "_".join([tss, args_tss.program]) + ".gff"
            predict = os.path.join(self.gff_outfolder, filename)
            print("Running merge and classify manual ....")
            stat_file = "stat_compare_TSSpredator_manual_{0}.csv".format(tss)
            merge_manual_predict_tss(
                predict, stat_file,
                os.path.join(self.tmps["tss"], filename),
                os.path.join(args_tss.gffs, gff), args_tss)
            shutil.move(stat_file, os.path.join(args_tss.out_folder,
                                                "statistics", tss, stat_file))
        self.helper.move_all_content(self.tmps["tss"],
                                     self.gff_outfolder, [".gff"])
        shutil.rmtree(self.tmps["tss"])

    def _validate(self, tsss, args_tss):
        '''validate TSS with genome annotation'''
        print("Running validation of annotation....")
        for tss in tsss:
            for gff in os.listdir(args_tss.gffs):
                if (gff[:-4] == tss) and (".gff" in gff):
                    break
            stat_file = os.path.join(
                    self.stat_outfolder, tss,
                    "".join(["stat_gene_vali_", tss, ".csv"]))
            out_cds_file = os.path.join(args_tss.out_folder, "tmp.gff")
            if args_tss.program.lower() == "tss":
                compare_file = os.path.join(self.gff_outfolder,
                                            "_".join([tss, "TSS.gff"]))
            elif args_tss.program.lower() == "processing":
                compare_file = os.path.join(self.gff_outfolder,
                                            "_".join([tss, "processing.gff"]))
            validate_gff(compare_file, os.path.join(args_tss.gffs, gff),
                         stat_file, out_cds_file, args_tss.utr_length,
                         args_tss.program.lower())
            shutil.move(out_cds_file, os.path.join(args_tss.gffs, gff))

    def _compare_ta(self, tsss, args_tss):
        '''compare TSS with transcript'''
        detect = False
        print("Running compare transcript assembly and TSS ...")
        self.multiparser.parser_gff(args_tss.ta_files, "transcript")
        self.multiparser.combine_gff(args_tss.gffs, self.tmps["ta"],
                                     None, "transcript")
        for tss in tsss:
            stat_out = os.path.join(
                    self.stat_outfolder, tss, "".join([
                        "stat_compare_TSS_transcript_",
                        tss, ".csv"]))
            for ta in os.listdir(self.tmps["ta"]):
                filename = ta.split("_transcript")
                if (filename[0] == tss) and (filename[1] == ".gff"):
                    detect = True
                    break
            compare_file = os.path.join(self.gff_outfolder,
                                        "_".join([tss, "TSS.gff"]))
            if detect:
                stat_ta_tss(os.path.join(self.tmps["ta"], ta), compare_file,
                            stat_out, self.tmps["ta_tss"],
                            self.tmps["tss_ta"], args_tss.fuzzy)
                self.helper.sort_gff(self.tmps["tss_ta"], compare_file)
                self.helper.sort_gff(self.tmps["ta_tss"],
                                     os.path.join(args_tss.ta_files, ta))
                os.remove(self.tmps["tss_ta"])
                os.remove(self.tmps["ta_tss"])
                detect = False

    def _stat_tss(self, tsss, feature):
        print("Running statistaics.....")
        for tss in tsss:
            compare_file = os.path.join(self.gff_outfolder,
                                        "_".join([tss, feature]) + ".gff")
            stat_tsspredator(
                compare_file, feature,
                os.path.join(self.stat_outfolder, tss, "_".join([
                    "stat", feature, "class", tss]) + ".csv"),
                os.path.join(self.stat_outfolder, tss, "_".join([
                    "stat", feature, "libs", tss]) + ".csv"))
            self.helper.move_all_content(os.getcwd(), os.path.join(
                self.stat_outfolder, tss), ["_class", ".png"])
            if os.path.exists(os.path.join(
                    self.stat_outfolder, "TSSstatistics.tsv")):
                shutil.move(
                    os.path.join(
                        self.stat_outfolder, "TSSstatistics.tsv"),
                    os.path.join(
                        self.stat_outfolder, tss, "TSSstatistics.tsv"))
            plot_venn(compare_file, feature)
            self.helper.move_all_content(os.getcwd(), os.path.join(
                self.stat_outfolder, tss), ["_venn", ".png"])

    def _set_gen_config(self, args_tss, input_folder):
        prefixs = []
        detect = False
        for fasta in os.listdir(self.fasta_path):
            for gff in os.listdir(self.gff_path):
                if fasta[:-3] == gff[:-4]:
                    prefix = fasta[:-3]
                    for wig in os.listdir(self.wig_path):
                        filename = wig.split("_STRAIN_")
                        if filename[1][:-4] == prefix:
                            detect = True
                            break
                    if detect:
                        prefixs.append(prefix)
                        config = os.path.join(
                                input_folder,
                                "_".join(["config", prefix]) + ".ini")
                        self._gen_config(
                            prefix, args_tss,
                            os.path.join(self.gff_path, gff), self.wig_path,
                            os.path.join(self.fasta_path, fasta), config)
        return prefixs

    def _merge_wigs(self, wig_folder, prefix, libs):
        self.helper.check_make_folder(os.path.join(os.getcwd(),
                                      self.tmps["tmp"]))
        for wig_file in os.listdir(wig_folder):
            for lib in libs:
                info = lib.split(":")
                if (info[0][:-4] in wig_file) and (info[-1] == "+") and (
                        prefix in wig_file) and (
                        os.path.isfile(os.path.join(wig_folder, wig_file))):
                    Helper().merge_file(
                            os.path.join(wig_folder, wig_file),
                            os.path.join("tmp", "merge_forward.wig"))
                if (info[0][:-4] in wig_file) and (info[-1] == "-") and (
                        prefix in wig_file) and (
                        os.path.isfile(os.path.join(wig_folder, wig_file))):
                    Helper().merge_file(
                            os.path.join(wig_folder, wig_file),
                            os.path.join("tmp", "merge_reverse.wig"))

    def _check_orphan(self, prefixs, wig_folder, args_tss):
        '''if genome has no locus tag, it can use for classify the TSS'''
        for prefix in prefixs:
            self._merge_wigs(wig_folder, prefix, args_tss.libs)
            tmp_tss = os.path.join(self.tmps["tmp"], "_".join([
                          prefix, args_tss.program + ".gff"]))
            pre_tss = os.path.join(self.gff_outfolder, "_".join([
                          prefix, args_tss.program + ".gff"]))
            check_orphan(pre_tss, os.path.join(
                args_tss.gffs, prefix + ".gff"),
                "tmp/merge_forward.wig", "tmp/merge_reverse.wig", tmp_tss)
            shutil.move(tmp_tss, pre_tss)
        shutil.rmtree("tmp")

    def _remove_files(self, args_tss):
        print("Remove temperary files and folders...")
        self.helper.remove_tmp(args_tss.fastas)
        self.helper.remove_tmp(args_tss.gffs)
        self.helper.remove_tmp(args_tss.wig_folder)
        self.helper.remove_tmp(args_tss.ta_files)
        if "merge_forward.wig" in os.listdir(os.getcwd()):
            os.remove("merge_forward.wig")
        if "merge_reverse.wig" in os.listdir(os.getcwd()):
            os.remove("merge_reverse.wig")

    def _deal_with_overlap(self, out_folder, args_tss):
        '''deal with the situation that TSS and 
        processing site at the same position'''
        if args_tss.overlap_feature.lower() == "both":
            pass
        else:
            print("Comparing TSS and Processing site...")
            if args_tss.program.lower() == "tss":
                for tss in os.listdir(out_folder):
                    if tss.endswith("_TSS.gff"):
                        ref = self.helper.get_correct_file(
                                args_tss.references, "_processing.gff",
                                tss.replace("_TSS.gff", ""), None, None)
                        filter_tss_pro(os.path.join(out_folder, tss),
                                       ref, args_tss.overlap_feature,
                                       args_tss.cluster)
            elif args_tss.program.lower() == "processing_site":
                for tss in os.listdir(out_folder):
                    if tss.endswith("_processing.gff"):
                        ref = self.helper.get_correct_file(
                                args_tss.references, "_TSS.gff",
                                tss.replace("_processing.gff", ""), None, None)
                        filter_tss_pro(os.path.join(out_folder, tss),
                                       ref, args_tss.overlap_feature,
                                       args_tss.cluster)

    def _low_expression(self, args_tss, gff_folder):
        '''deal with the low expressed TSS'''
        prefix = None
        self._merge_wigs(args_tss.wig_folder, "wig", args_tss.libs)
        for gff in os.listdir(gff_folder):
            if (args_tss.program.lower() == "tss") and (
                    gff.endswith("_TSS.gff")):
                prefix = gff.replace("_TSS.gff", "")
            elif (args_tss.program.lower() == "processing") and (
                    gff.endswith("_processing.gff")):
                prefix = gff.replace("_processing.gff", "")
            if prefix:
                out = open(os.path.join(
                    self.stat_outfolder, prefix, "_".join([
                        "stat", prefix, "low_expression_cutoff.csv"])), "w")
                out.write("\t".join(["strain", "cutoff_coverage"]) + "\n")
                cutoff = filter_low_expression(
                        os.path.join(gff_folder, gff), args_tss,
                        "tmp/merge_forward.wig", "tmp/merge_reverse.wig",
                        "tmp/without_low_expression.gff")
                out.write("\t".join([prefix, str(cutoff)]) + "\n")
                os.remove(os.path.join(gff_folder, gff))
                shutil.move("tmp/without_low_expression.gff",
                            os.path.join(gff_folder, gff))
                prefix = None
        out.close()

    def run_tsspredator(self, args_tss):
        input_folder = os.path.join(args_tss.out_folder, "configs")
        for gff in os.listdir(args_tss.gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(
                                                 args_tss.gffs, gff))
        self.helper.check_make_folder(self.gff_outfolder)
        self.multiparser.parser_fasta(args_tss.fastas)
        self.multiparser.parser_gff(args_tss.gffs, None)
        self.multiparser.parser_wig(args_tss.wig_folder)
        prefixs = self._set_gen_config(args_tss, input_folder)
        for prefix in prefixs:
            out_path = os.path.join(
                    self.master, "_".join(["MasterTable", prefix]))
            config_file = os.path.join(
                    input_folder, "_".join(["config", prefix]) + ".ini")
            self._start_to_run(args_tss.tsspredator_path, config_file,
                               out_path, prefix)
            if os.path.exists(os.path.join(out_path, "TSSstatistics.tsv")):
                shutil.move(os.path.join(out_path, "TSSstatistics.tsv"),
                            os.path.join(
                                self.stat_outfolder, "TSSstatistics.tsv"))
        if args_tss.program.lower() == "processing_site":
            args_tss.program = "processing"
        self._convert_gff(prefixs, args_tss)
        if args_tss.check_orphan:
            print("checking the orphan TSS...")
            self._check_orphan(prefixs,
                               os.path.join(args_tss.wig_folder, "tmp"),
                               args_tss)
        self.multiparser.combine_gff(args_tss.gffs, self.gff_outfolder,
                                     None, args_tss.program)
        datas = []
        for gff in os.listdir(self.gff_outfolder):
            if gff.endswith(".gff"):
                gff_folder = gff.replace("".join(["_", args_tss.program,
                                                  ".gff"]), "")
                self.helper.check_make_folder(
                     os.path.join(self.stat_outfolder, gff_folder))
                datas.append(gff_folder)
        if args_tss.remove_low_expression is not None:
            self._low_expression(args_tss, self.gff_outfolder)
        if args_tss.manual is not None:
            self.multiparser.combine_wig(args_tss.gffs, self.wig_path,
                                         None, args_tss.libs)
            self._merge_manual(datas, args_tss)
        self._deal_with_overlap(self.gff_outfolder, args_tss)
        if args_tss.stat:
            self._stat_tss(datas, args_tss.program)
        if args_tss.validate:
            self._validate(datas, args_tss)
        if args_tss.ta_files is not None:
            self._compare_ta(datas, args_tss)
        self._remove_files(args_tss)
