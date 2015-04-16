#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
import shutil
from subprocess import call, Popen
from annogesiclib.helper import Helper
from annogesiclib.multiparser import Multiparser
from annogesiclib.converter import Converter
from annogesiclib.merge_manual import merge_manual_predict_tss
from annogesiclib.stat_TSSpredator import stat_tsspredator
from annogesiclib.plot_TSS_venn import plot_venn
from annogesiclib.validate_gene import validate_gff
from annogesiclib.stat_TA_comparison import stat_ta_tss

class TSSpredator(object):

    def __init__(self, out_folder, TA_files, gffs, wig_folder, fastas):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()
        self.master = os.path.join(out_folder, "MasterTables")
        self.tmp_tss = "tmp_TSS"
        self.tmp_ta = os.path.join(TA_files, "tmp")
        self.tmp_ta_tss = "tmp_ta_tss"
        self.tmp_tss_ta = "tmp_tss"
        self.gff_path = os.path.join(gffs, "tmp")
        self.wig_path = os.path.join(wig_folder, "tmp")
        self.fasta_path = os.path.join(fastas, "tmp")
        self.stat_outfolder = os.path.join(out_folder, "statistics")
        self.gff_outfolder = os.path.join(out_folder, "gffs")

    def _assign_dict(self, lib_datas):
        return {"wig": lib_datas[0],
                "tex": lib_datas[1],
                "condition": int(lib_datas[2]),
                "replicate": lib_datas[3],
                "strand": lib_datas[4]}

    def _print_lib(self, lib_num, lib_list, out, wig_folder, prefix):
        for num_id in range(1, lib_num+1):
            cond_list = []
            for lib in lib_list:
                if num_id == lib["condition"]:
                    cond_list.append(lib)
            cond_sort_list = sorted(cond_list, key=lambda k: k['replicate'])
            for cond in cond_sort_list:
                out.write("{0}_{1}{2} = {3}\n".format(
                          prefix, cond["condition"], cond["replicate"],
                          os.path.join(wig_folder, cond["wig"])))

    def _start_to_run(self, tsspredator_path, config_file, out_path, prefix):
        print("Running TSSpredator for " + prefix)
        out = open(os.path.join(out_path, "log.txt"), "w")
        err = open(os.path.join(out_path, "err.txt"), "w")
        call(["java",  "-jar", tsspredator_path, 
              config_file], stdout=out ,stderr=err)
        out.close()
        err.close()

    def _import_lib(self, libs, wig_folder, project_strain_name, rep_set, lib_dict,
                    out, gff, program, list_num_id, fasta):
        lib_num = 0
        print("Runniun {0} now...".format(program))
        for lib in libs:
            lib_datas = lib.split(":")
            if lib_datas[0].endswith(".wig") is not True:
                print("Error:Exist a not proper wig files!!")
                sys.exit()
            for wig in os.listdir(wig_folder):
                filename = wig.split("_STRAIN_")
                if (filename[0] == lib_datas[0][:-4]) and \
                   (filename[1][:-4] == project_strain_name):
                    lib_datas[0] = wig
            if int(lib_datas[2]) > lib_num:
                lib_num = int(lib_datas[2])
            if lib_datas[3] not in rep_set:
                rep_set.add(lib_datas[3])
            if (lib_datas[1] == "tex") and \
               (lib_datas[4] == "+"):
                lib_dict["fp"].append(self._assign_dict(lib_datas))
            elif (lib_datas[1] == "tex") and \
                 (lib_datas[4] == "-"):
                lib_dict["fm"].append(self._assign_dict(lib_datas))
            elif (lib_datas[1] == "notex") and \
                 (lib_datas[4] == "+"):
                lib_dict["np"].append(self._assign_dict(lib_datas))
            elif (lib_datas[1] == "notex") and \
                 (lib_datas[4] == "-"):
                lib_dict["nm"].append(self._assign_dict(lib_datas))
        for num_id in range(1, lib_num+1):
            out.write("annotation_{0} = {1}\n".format(num_id, gff))
        if program == "TSS":
            self._print_lib(lib_num, lib_dict["fm"], out, wig_folder, "fivePrimeMinus")
            self._print_lib(lib_num, lib_dict["fp"], out, wig_folder, "fivePrimePlus")
        elif program == "processing_site":
            self._print_lib(lib_num, lib_dict["nm"], out, wig_folder, "fivePrimeMinus")
            self._print_lib(lib_num, lib_dict["np"], out, wig_folder, "fivePrimePlus")
        else:
            print("Error: Wrong program name!!!")
            sys.exit()
        for num_id in range(1, lib_num+1):
            out.write("genome_{0} = {1}\n".format(num_id, fasta))
        for num_id in range(1, lib_num+1):
            list_num_id.append(str(num_id))
        return (lib_num, num_id)

    def _gen_config(self, project_strain_name, out_folder,
                    libs, gff, wig_folder, fasta, height,
                    factor, height_reduction, factor_reduction,
                    base_height, output_prefixs, config_file, 
                    program, repmatch, cluster, utr_length):
        lib_dict = {"fp": [], "fm": [], "nm": [], "np": []}
        rep_set = set()
        master_folder = "MasterTable_" + project_strain_name
        out_path = os.path.join(self.master, master_folder)
        self.helper.check_make_folder(out_path)
        out = open(config_file, "w")
        out.write("TSSinClusterSelectionMethod = HIGHEST\n")
        out.write("allowedCompareShift = 1\n")
        out.write("allowedRepCompareShift = 1\n")
        list_num_id = []
        nums = self._import_lib(libs, wig_folder, project_strain_name, rep_set, lib_dict,
                                out, gff, program, list_num_id, fasta)
        lib_num = nums[0]
        num_id = nums[1]
        out.write("idList = ")
        out.write(",".join(list_num_id) + "\n")
        out.write("maxASutrLength = 100\n")
        out.write("maxGapLengthInGene = 500\n")
        out.write("maxNormalTo5primeFactor = 1.5\n")
        out.write("maxTSSinClusterDistance = {0}\n".format(cluster + 1))
        out.write("maxUTRlength = {0}\n".format(utr_length))
        out.write("min5primeToNormalFactor = 2\n")
        out.write("minCliffFactor = {0}\n".format(factor))
        out.write("minCliffFactorDiscount = {0}\n".format(factor_reduction))
        out.write("minCliffHeight = {0}\n".format(height))
        out.write("minCliffHeightDiscount = {0}\n".format(height_reduction))
        out.write("minNormalHeight = {0}\n".format(base_height))
        out.write("minNumRepMatches = {0}\n".format(repmatch))
        out.write("minPlateauLength = 0\n")
        out.write("mode = cond\n")
        out.write("normPercentile = 0.9\n")
        if program == "TSS":
            self._print_lib(lib_num, lib_dict["nm"], out, wig_folder,"normalMinus")
            self._print_lib(lib_num, lib_dict["np"], out, wig_folder,"normalPlus")
        else:
            self._print_lib(lib_num, lib_dict["fm"], out, wig_folder,"normalMinus")
            self._print_lib(lib_num, lib_dict["fp"], out, wig_folder,"normalPlus")
        out.write("numReplicates = {0}\n".format(len(rep_set)))
        out.write("numberOfDatasets = {0}\n".format(lib_num))
        out.write("outputDirectory = {0}\n".format(out_path))
        for prefix_id in range(len(output_prefixs)):
            out.write("outputPrefix_{0} = {1}\n".format(prefix_id + 1, output_prefixs[prefix_id]))
        out.write("projectName = {0}\n".format(project_strain_name))
        out.write("superGraphCompatibility = igb\n")
        out.write("texNormPercentile = 0.5\n")
        out.write("writeGraphs = 0\n")
        out.write("writeNocornacFiles = 0\n")
        out.close()

    def _convert_gff(self, prefixs, out_folder, feature):
        for prefix in prefixs:
            out_file = os.path.join(self.gff_outfolder, "_".join([prefix, feature]) + ".gff")
            gff_f = open(out_file, "w")
            out_path = os.path.join(self.master, "_".join(["MasterTable", prefix]))
            if "MasterTable.tsv" not in os.listdir(out_path):
                print("Error:there is not MasterTable file in {0}".format(out_path))
                print("Please check config.ini")
            else:
                self.converter.convert_mastertable2gff(
                                os.path.join(out_path, "MasterTable.tsv"),
                                "TSSpredator", feature, prefix, out_file)
            gff_f.close()

    def _merge_manual(self, tsss, gffs, manual, wig_path,
                      out_folder, cluster, feature, length, libs):
        self.helper.check_make_folder(os.path.join(os.getcwd(), self.tmp_tss))
        for tss in tsss:
            for gff in os.listdir(gffs):
                if (gff[:-4] == tss) and (".gff" in gff):
                    break
            filename = "_".join([tss, feature]) + ".gff"
            predict = os.path.join(self.gff_outfolder, filename)
            print("Running merge and classify manual ....")
            stat_file = "stat_compare_TSSpredator_manual_{0}.csv".format(tss)
            merge_manual_predict_tss(predict, manual, stat_file, 
                            os.path.join(self.tmp_tss, filename),
                            os.path.join(gffs, gff), cluster, length, libs,
                            wig_path, feature)
            os.rename(stat_file, os.path.join(out_folder, "statistics", tss, stat_file))
        self.helper.move_all_content(self.tmp_tss, self.gff_outfolder, [".gff"])
        shutil.rmtree(self.tmp_tss)

    def _validate(self, tsss, gffs, utr_length, out_folder):
        print("Running validation of annotation....")
        for tss in tsss:
            for gff in os.listdir(gffs):
                if (gff[:-4] == tss) and (".gff" in gff):
                    break
            stat_file = os.path.join(self.stat_outfolder, tss, 
                        "".join(["stat_gene_vali_", tss, ".csv"]))
            out_cds_file = os.path.join(out_folder, "tmp.gff")
            compare_file = os.path.join(self.gff_outfolder, "_".join([tss, "TSS.gff"]))
            validate_gff(compare_file, os.path.join(gffs, gff), 
                         stat_file, out_cds_file, utr_length)
            os.rename(out_cds_file, os.path.join(gffs, gff))

    def _compare_ta(self, TA_files, gffs, tsss, fuzzy):
        detect = False
        print("Running compare transcript assembly and TSS ...")
        self.multiparser._parser_gff(TA_files, "transcript")
        self.multiparser._combine_gff(gffs, self.tmp_ta, None, "transcript")
        for tss in tsss:
            stat_out = os.path.join(self.stat_outfolder, tss, 
                       "".join(["stat_compare_TSS_Transcriptome_assembly_", tss, ".csv"]))
            for ta in os.listdir(self.tmp_ta):
                filename = ta.split("_transcript")
                if (filename[0] == tss) and (filename[1] == ".gff"):
                    detect = True
                    break
            compare_file = os.path.join(self.gff_outfolder, "_".join([tss, "TSS.gff"]))
            if detect:
                stat_ta_tss(os.path.join(self.tmp_ta, ta), compare_file, 
                            stat_out, self.tmp_ta_tss, self.tmp_tss_ta, fuzzy)
                self.helper.sort_gff(self.tmp_tss_ta, compare_file)
                self.helper.sort_gff(self.tmp_ta_tss, os.path.join(TA_files, ta))
                os.remove(self.tmp_tss_ta)
                os.remove(self.tmp_ta_tss)
                detect = False

    def _stat_tss(self, tsss, manual, feature):
        print("Running statistaics.....")
        for tss in tsss:
            compare_file = os.path.join(self.gff_outfolder, 
                           "_".join([tss, feature]) + ".gff")
            stat_tsspredator(compare_file, feature, 
                             os.path.join(self.stat_outfolder, tss, 
                             "_".join(["stat", feature, "class", tss]) + ".csv"))
            self.helper.move_all_content(os.getcwd(), 
                        os.path.join(self.stat_outfolder, tss), ["_class", ".png"])
            #### generate venn diagram
            plot_venn(compare_file, feature)
            self.helper.move_all_content(os.getcwd(),
                        os.path.join(self.stat_outfolder, tss), ["_venn", ".png"])

    def _set_gen_config(self, feature, input_folder, out_folder, libs, height, factor,
                        height_reduction, factor_reduction, utr_length, base_height, 
                        output_prefixs, repmatch, cluster):
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
                    if detect is True:
                        prefixs.append(prefix)
                        config = os.path.join(input_folder, "_".join(["config", prefix]) + ".ini")
                        self._gen_config(prefix, out_folder,
                                         libs, os.path.join(self.gff_path, gff), self.wig_path,
                                         os.path.join(self.fasta_path, fasta), height,
                                         factor, height_reduction, factor_reduction,
                                         base_height, output_prefixs, config, feature, 
                                         repmatch, cluster, utr_length)
        return prefixs

    def run_tsspredator(self, tsspredator_path, program, input_folder, fastas, gffs, 
                        wig_folder, libs, output_prefixs, height, height_reduction, 
                        factor, factor_reduction, base_height, repmatch, out_folder, 
                        project_path, stat, validate, manual, TA_files, fuzzy, 
                        utr_length, cluster, nt_length):
        ####  First of all, generate config file for running.
        for gff in os.listdir(gffs):
            if gff.endswith(".gff"):
                self.helper.check_uni_attributes(os.path.join(gffs, gff))
        self.helper.check_make_folder(self.gff_outfolder)
        self.multiparser._parser_fasta(fastas)
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._parser_wig(wig_folder)
        prefixs = self._set_gen_config(program, input_folder, out_folder, libs, height, 
                                       factor, height_reduction, factor_reduction, utr_length, 
                                       base_height, output_prefixs, repmatch, cluster)
        #### After generating config file, we can start to run TSSpredator.
        for prefix in prefixs:
            out_path = os.path.join(self.master, "_".join(["MasterTable", prefix]))
            config_file = os.path.join(input_folder, "_".join(["config", prefix]) + ".ini")
            self._start_to_run(tsspredator_path, config_file, out_path, prefix)
        if program == "processing_site":
            program = "processing"
        self._convert_gff(prefixs, out_folder, program)
        #### based on the original file to merge the information of strains.
        self.multiparser._combine_gff(gffs, self.gff_outfolder, None, program)
        datas = []
        for gff in os.listdir(self.gff_outfolder):
            if gff.endswith(".gff"):
                gff_folder = gff.replace("".join(["_", program, ".gff"]), "")
                self.helper.check_make_folder(os.path.join(self.stat_outfolder, gff_folder))
                datas.append(gff_folder)
        if manual is not False:
            self.multiparser._combine_wig(gffs, self.wig_path, None)
            self._merge_manual(datas, gffs, manual, wig_folder, out_folder, 
                               cluster, program, nt_length, libs)
        if stat is not False:
             self._stat_tss(datas, manual, program)
        if validate is not False:
            self._validate(datas, gffs, utr_length, out_folder)
        if TA_files is not False:
            self._compare_ta(TA_files, gffs, datas, fuzzy)
        #### Remove temperary folders and files.
        print("Remove temperary files and folders...")
        self.helper.remove_tmp(fastas)
        self.helper.remove_tmp(gffs)
        self.helper.remove_tmp(wig_folder)
        self.helper.remove_tmp(TA_files)
