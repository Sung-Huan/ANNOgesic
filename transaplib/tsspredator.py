#!/usr/bin/python
from Bio import SeqIO
import os	
import sys
from subprocess import call, Popen
from transaplib.helper import Helper
from transaplib.multiparser import Multiparser
from transaplib.converter import Converter
from transaplib.merge_manual import Merge_manual_predict_TSS 
from transaplib.stat_TSSpredator import Stat_TSSpredator
from transaplib.plot_TSS_venn import Plot_Venn
from transaplib.validate_gene import Validate_gff
from transaplib.stat_TA_comparison import Stat_TA_TSS

class TSSpredator(object):

    def __init__(self):
        self.multiparser = Multiparser()
        self.helper = Helper()
        self.converter = Converter()

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
                out.write("%s_%s%s = %s%s\n" % 
                         (prefix, str(cond["condition"]), cond["replicate"],
                          wig_folder, cond["wig"]))

    def _start_to_run(self, tsspredator_path, config_file, out_path, prefix):
        print("Running TSSpredator for " + prefix)
        out = open(out_path + "/log.txt", "w")
        err = open(out_path + "/err.txt", "w")
        call(["java",  "-jar", tsspredator_path, 
              config_file], stdout=out ,stderr=err)
        out.close()
        err.close()

    def _import_lib(self, libs, wig_folder, project_strain_name, rep_set, lib_dict,
                    out, gff, program, list_num_id, fasta):
        lib_num = 0
        print(program)
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
            out.write("annotation_%s = %s\n" % (str(num_id), gff))
        if program == "TSS":
            self._print_lib(lib_num, lib_dict["fm"], out, wig_folder, "fivePrimeMinus")
            self._print_lib(lib_num, lib_dict["fp"], out, wig_folder,"fivePrimePlus")
        elif program == "processing_site":
            self._print_lib(lib_num, lib_dict["nm"], out, wig_folder, "fivePrimeMinus")
            self._print_lib(lib_num, lib_dict["np"], out, wig_folder,"fivePrimePlus")
        else:
            print("Error: Wrong program name!!!")
            sys.exit()
        for num_id in range(1, lib_num+1):
            out.write("genome_%s = %s\n" % (str(num_id), fasta))
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
        out_path =  out_folder + "/MasterTables/MasterTable_" + project_strain_name
        if "MasterTable_" + project_strain_name in os.listdir(out_folder + "/MasterTables"):
            call(["rm", "-rf", out_path])
        call(["mkdir", out_path])
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
        out.write("maxTSSinClusterDistance = %s\n" % str(cluster + 1))
        out.write("maxUTRlength = {}\n".format(utr_length))
        out.write("min5primeToNormalFactor = 2\n")
        out.write("minCliffFactor = %s\n" % factor)
        out.write("minCliffFactorDiscount = %s\n" % factor_reduction)
        out.write("minCliffHeight = %s\n" %height)
        out.write("minCliffHeightDiscount = %s\n" % height_reduction)
        out.write("minNormalHeight = %s\n" % base_height)
        out.write("minNumRepMatches = %s\n" % repmatch)
        out.write("minPlateauLength = 0\n")
        out.write("mode = cond\n")
        out.write("normPercentile = 0.9\n")
        if program == "TSS":
            self._print_lib(lib_num, lib_dict["nm"], out, wig_folder,"normalMinus")
            self._print_lib(lib_num, lib_dict["np"], out, wig_folder,"normalPlus")
        else:
            self._print_lib(lib_num, lib_dict["fm"], out, wig_folder,"normalMinus")
            self._print_lib(lib_num, lib_dict["fp"], out, wig_folder,"normalPlus")
        out.write("numReplicates = %s\n" % str(len(rep_set)))
        out.write("numberOfDatasets = %s\n" % lib_num)
        out.write("outputDirectory = %s\n" % out_path)
        for prefix_id in range(len(output_prefixs)):
            out.write("outputPrefix_%s = %s\n" % (str(prefix_id + 1), output_prefixs[prefix_id]))
        out.write("projectName = %s\n" % project_strain_name)
        out.write("superGraphCompatibility = igb\n")
        out.write("texNormPercentile = 0.5\n")
        out.write("writeGraphs = 0\n")
        out.write("writeNocornacFiles = 0\n")
        out.close()

    def _convert_gff(self, prefixs, gff_outfolder, out_folder, feature):
        for prefix in prefixs:
            gff_f = open(gff_outfolder + prefix + "_" + feature + ".gff", "w")
            out_path = out_folder + "/MasterTables/MasterTable_" + prefix
            if "MasterTable.tsv" not in os.listdir(out_path):
                print("Error:there is not MasterTable file in " + out_path)
                print("Please check config.ini")
            else:
                self.converter.Convert_Mastertable2gff(
                                out_path + "/MasterTable.tsv",
                                "TSSpredator", feature, prefix,
                                gff_outfolder + prefix + "_" + feature + ".gff")
            gff_f.close()

    def _merge_manual(self, tsss, gffs, gff_outfolder, manual, wig_path,
                      out_folder, cluster, feature, length, libs):
        call(["mkdir", "tmp_TSS"])
        for tss in tsss:
            for gff in os.listdir(gffs):
                if (gff[:-4] == tss) and (".gff" in gff):
                    break
            predict = gff_outfolder + tss + "_" + feature + ".gff"
            print("Running merge and classify manual ....")
            stat_file = "stat_compare_TSSpredator_manual_" + tss + ".csv"
            Merge_manual_predict_TSS(predict, manual, stat_file, 
                                     "tmp_TSS/" + tss + "_" + feature + ".gff",
                                     gffs + "/" + gff, cluster, length, libs,
                                     wig_path, feature)
            call(["mv", stat_file, out_folder + "/statistics/" + tss])
        os.system("mv tmp_TSS/*.gff " + gff_outfolder)
        call(["rm", "-rf", "tmp_TSS"])

    def _validate(self, tsss, gffs, utr_length, out_folder, gff_outfolder, stat_outfolder):
        print("Running validation of annotation....")
        for tss in tsss:
            for gff in os.listdir(gffs):
                if (gff[:-4] == tss) and (".gff" in gff):
                    break
            stat_file = stat_outfolder + tss + "/stat_gene_vali_" + tss + ".csv"
            out_cds_file = out_folder + "/tmp.gff"
            compare_file = gff_outfolder + tss + "_TSS.gff"
            Validate_gff(compare_file, gffs + "/" + gff, stat_file, out_cds_file, utr_length)
            call(["mv", out_cds_file, gffs + "/" + gff])

    def _compare_TA(self, TA_files, gffs, tsss, stat_outfolder, gff_outfolder, fuzzy):
        detect = False
        print("Running compare transcript assembly and TSS ...")
        self.multiparser._parser_gff(TA_files, "transcript")
        self.multiparser._combine_gff(gffs, TA_files + "/tmp", None, "transcript")
        for tss in tsss:
            stat_out = stat_outfolder + tss + "/stat_compare_TSS_Transcriptome_assembly_" + tss + ".csv"
            for ta in os.listdir(TA_files + "/tmp"):
                filename = ta.split("_transcript")
                if (filename[0] == tss) and (filename[1] == ".gff"):
                    detect = True
                    break
            compare_file = gff_outfolder + tss + "_TSS.gff"
            if detect:
                Stat_TA_TSS(TA_files + "/tmp/" + ta, compare_file, 
                            stat_out, "tmp_ta_tss", "tmp_tss", fuzzy)
                self.helper.sort_gff("tmp_tss", compare_file)
                self.helper.sort_gff("tmp_ta_tss", TA_files + "/" + ta)
                call(["rm", "tmp_tss"])
                call(["rm", "tmp_ta_tss"])
                detect = False

    def _stat_TSS(self, tsss, manual, gff_outfolder, 
                  stat_outfolder, feature):
        print("Running statistaics.....")
        for tss in tsss:
            compare_file = gff_outfolder + tss + "_" + feature + ".gff"
            Stat_TSSpredator(compare_file, feature, 
                             "".join([stat_outfolder + tss + \
                                      "/stat_" + feature + \
                                      "_class_" + tss + ".csv"]))
            os.system("mv " + feature + "_class*.png " + stat_outfolder + tss)
            #### generate venn diagram
            Plot_Venn(compare_file, feature)
            os.system("mv " + feature + "_venn*.png " + stat_outfolder + tss)

    def _set_gen_config(self, fasta_path, gff_path, feature, input_folder,
                        out_folder, libs, wig_path, height, factor,
                        height_reduction, factor_reduction, utr_length, base_height, 
                        output_prefixs, repmatch, cluster):
        prefixs = []
        detect = False
        for fasta in os.listdir(fasta_path):
            for gff in os.listdir(gff_path):
                if fasta[:-3] == gff[:-4]:
                    prefix = fasta[:-3]
                    for wig in os.listdir(wig_path):
                        filename = wig.split("_STRAIN_")
                        if filename[1][:-4] == prefix:
                            detect = True
                            break
                    if detect is True:
                        prefixs.append(prefix)
                        config = input_folder + "/config_" + prefix + ".ini"
                        self._gen_config(prefix, out_folder,
                                         libs, gff_path + gff, wig_path,
                                         fasta_path + fasta, height,
                                         factor, height_reduction, factor_reduction,
                                         base_height, output_prefixs, config, feature, 
                                         repmatch, cluster, utr_length)
        return prefixs

    def run_TSSpredator(self, tsspredator_path, input_folder, fastas, gffs, 
                        wig_folder, libs, output_prefixs, height, height_reduction, 
                        factor, factor_reduction, base_height, repmatch, out_folder, 
                        project_path, stat, validate, manual, TA_files, fuzzy, 
                        utr_length, cluster, nt_length):
        ####  First of all, generate config file for running.
#        for gff in os.listdir(gffs):
#            if gff.endswith(".gff"):
#                self.helper.check_uni_attributes(gffs + "/" + gff)
        gff_outfolder = out_folder + "/gffs/"
        os.system("rm -rf " + gff_outfolder + "*")
        self.multiparser._parser_fasta(fastas)
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._parser_wig(wig_folder)
        gff_path = gffs + "/tmp/"
        wig_path = wig_folder + "/tmp/"
        fasta_path = fastas + "/tmp/"
        stat_outfolder = out_folder + "/statistics/"
        merge_wigs = os.getcwd() + "/merge_wigs"
        prefixs = self._set_gen_config(fasta_path, gff_path, "TSS", input_folder, 
                                       out_folder, libs, wig_path, height, factor,
                                       height_reduction, factor_reduction, utr_length, 
                                       base_height, output_prefixs, repmatch, cluster)
        #### After generating config file, we can start to run TSSpredator.
        for prefix in prefixs:
            out_path = out_folder + "/MasterTables/MasterTable_" + prefix
            config_file = input_folder + "/config_" + prefix + ".ini"
            self._start_to_run(tsspredator_path, config_file, out_path, prefix)
        self._convert_gff(prefixs, gff_outfolder, out_folder, "TSS")
        #### based on the original file to merge the information of strains.
        self.multiparser._combine_gff(gffs, gff_outfolder, None, "TSS")
        tsss = []
        for tss in os.listdir(gff_outfolder):
            if tss.endswith(".gff"):
                self.helper.check_make_folder(stat_outfolder, tss.replace("_TSS.gff", ""))
                tsss.append(tss.replace("_TSS.gff", ""))
        if manual is not False:
            self.multiparser._combine_wig(gffs, wig_path, None)
            self._merge_manual(tsss, gffs, gff_outfolder, manual, wig_folder,
                               out_folder, cluster, "TSS", nt_length, libs)
        if stat is not False:
             self._stat_TSS(tsss, manual, gff_outfolder, stat_outfolder, "TSS")
        if validate is not False:
            self._validate(tsss, gffs, utr_length, out_folder, gff_outfolder, stat_outfolder)
        if TA_files is not False:
            self._compare_TA(TA_files, gffs, tsss, stat_outfolder, gff_outfolder, fuzzy)
        #### Remove temperary folders and files.
#        print("Remove temperary files and folders...")
#        self.helper.remove_tmp(fastas)
#        self.helper.remove_tmp(gffs)
#        self.helper.remove_tmp(wig_folder)
#        self.helper.remove_tmp(TA_files)
    def run_processing_site(self, tsspredator_path, input_folder, fastas, 
                            gffs, wig_folder, libs, output_prefixs, height, 
                            height_reduction, factor, factor_reduction, 
                            base_height, repmatch, out_folder, project_path, 
                            stat, manual, cluster, utr_length, nt_length):
        #### First of all, generate config file for running.
#        for gff in os.listdir(gffs):
#            if gff.endswith(".gff"):
#                self.helper.check_uni_attributes(gffs + "/" + gff)
        self.multiparser._parser_fasta(fastas)
        self.multiparser._parser_gff(gffs, None)
        self.multiparser._parser_wig(wig_folder)
        gff_path = gffs + "/tmp/"
        wig_path = wig_folder + "/tmp/"
        fasta_path = fastas + "/tmp/"
        gff_outfolder = out_folder + "/gffs/"
        stat_outfolder = out_folder + "/statistics/"
        prefixs = self._set_gen_config(fasta_path, gff_path, 
                                       "processing_site", input_folder,
                                       out_folder, libs, wig_path, 
                                       height, factor, height_reduction, 
                                       factor_reduction, utr_length, base_height,
                                       output_prefixs, repmatch, cluster)
        #### After generating config file, we can start to run TSSpredator.
        for prefix in prefixs:
            out_path = out_folder + "/MasterTables/MasterTable_" + prefix
            config_file = input_folder + "/config_" + prefix + ".ini"
            self._start_to_run(tsspredator_path, config_file, out_path, prefix)
        self._convert_gff(prefixs, gff_outfolder, out_folder, "processing")
        #### based on the original file to merge the information of strains.
        self.multiparser._combine_gff(gffs, gff_outfolder, None, "processing")
        processings = []
        for processing in os.listdir(gff_outfolder):
            if ".gff" in processing:
                self.helper.check_make_folder(stat_outfolder, processing.replace("_processing.gff", ""))
                processings.append(processing.replace("_processing.gff", ""))
        if manual is not False:
            self._merge_manual(processings, gffs, gff_outfolder, manual, wig_folder,
                               out_folder, cluster, "processing", nt_length, libs)
        if stat is not False:
            self._stat_TSS(processings, manual, gff_outfolder, stat_outfolder, "processing")
        #### Remove temperary folders and files.
#        print("Remove temperary files and folders...")
#        self.helper.remove_tmp(fastas)
#        self.helper.remove_tmp(gffs)
#        self.helper.remove_tmp(wig_folder)
