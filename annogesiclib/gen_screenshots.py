import os
import sys
from glob import glob
from annogesiclib.gff3 import Gff3Parser
from annogesiclib.parser_wig import WigParser


def load_wigs(out, lib_t, lib_n, lib_f):
    if lib_t and lib_n:
        for index in range(len(lib_t)):
            out.write("load {0}\n".format(
                      os.path.join(os.getcwd(), lib_t[index])))
            out.write("load {0}\n".format(
                      os.path.join(os.getcwd(), lib_n[index])))
    elif lib_t:
        for lib in lib_t:
            out.write("load {0}\n".format(os.path.join(os.getcwd(), lib)))
    elif lib_n:
        for lib in lib_n:
            out.write("load {0}\n".format(os.path.join(os.getcwd(), lib)))
    if lib_f:
        for lib in lib_f:
            out.write("load {0}\n".format(os.path.join(os.getcwd(), lib)))


def set_data_range(out, gff, wigs, strand):
    '''set and print the DataRange'''
    max_range = 0
    for strains in wigs.values():
        for strain, wig_datas in strains.items():
            if strain == gff.seq_id:
                for wig in wig_datas[(gff.start - 1): gff.end]:
                    if max_range < wig.coverage:
                        max_range = wig.coverage
    max_range = int(max_range / 1) + 10
    if strand == "+":
        out.write("setDataRange 0,{0}\n".format(max_range))
    else:
        max_range = max_range * -1
        out.write("setDataRange 0,{0}\n".format(max_range))


def print_batch(args_sc, out, strand, lib_t, lib_n, lib_f, strain):
    '''print the batch file'''
    out.write("new\n")
    out.write("genome {0}\n".format(os.path.join(os.getcwd(), args_sc.fasta)))
    out.write("load {0}\n".format(os.path.join(os.getcwd(), args_sc.main_gff)))
    gff = args_sc.main_gff.split("/")
    out.write("{0} {1}\n".format(args_sc.present, gff[-1]))
    if args_sc.side_gffs is not None:
        for files in args_sc.side_gffs:
            for filename in glob(files):
                out.write("load {0}\n".format(os.path.join(os.getcwd(), filename)))
                gff = filename.split("/")
                out.write("{0} {1}\n".format(args_sc.present, gff[-1]))
    load_wigs(out, lib_t, lib_n, lib_f)
    out.write("maxPanelHeight {0}\n".format(args_sc.height))
    if strand == "+":
        out.write("snapshotDirectory {0}\n".format(
                      os.path.join(os.getcwd(), args_sc.output_folder,
                                   strain, "forward")))
    else:
        out.write("snapshotDirectory {0}\n".format(
                      os.path.join(os.getcwd(), args_sc.output_folder,
                                   strain, "reverse")))


def import_wig(lib, wigs, strand):
    wig_parser = WigParser()
    for wig in lib:
        wigs[wig] = {}
        strain = ""
        wig_fh = open(wig)
        for entry in wig_parser.parser(wig_fh, strand):
            if strain != entry.strain:
                wigs[wig][entry.strain] = []
                strain = entry.strain
            wigs[wig][strain].append(entry)
        wig_fh.close()


def gen_batch(lib_t, lib_n, lib_f, strand, gffs, out, seq):
    '''generate the batch file'''
    wigs = {}
    if lib_t and lib_n:
        import_wig(lib_t, wigs, strand)
        import_wig(lib_n, wigs, strand)
    elif lib_t:
        import_wig(lib_t, wigs, strand)
    elif lib_n:
        import_wig(lib_n, wigs, strand)
    if lib_f:
        import_wig(lib_f, wigs, strand)
    if strand == "+":
        print("Printing the forward batch files...")
    else:
        print("Printing the reverse batch files...")
    for gff in gffs:
        if gff.seq_id not in seq.keys():
            print("Error: The genome names in fasta file "
                  "and gff file are different!!")
            sys.exit()
        if (gff.start - 200) <= 0:
            start = 1
        else:
            start = gff.start - 200
        if (gff.end + 200) >= len(seq[gff.seq_id]):
            end = len(seq[gff.seq_id])
        else:
            end = gff.end + 200
        out.write("goto {0}:{1}-{2}\n".format(
                  gff.seq_id, start, end))
        set_data_range(out, gff, wigs, strand)
        out.write("snapshot {0}:{1}-{2}.png\n".format(
                  gff.seq_id, gff.start, gff.end))


def get_length(fasta_file):
    '''get sequence information and we can know the length of seq'''
    seq = {}
    with open(fasta_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                strain = line[1:]
                seq[strain] = ""
            else:
                seq[strain] = seq[strain] + line
    return seq


def gen_screenshot(args_sc, libs, forward_file, reverse_file, strain):
    '''Generation of screenshot of IGV for reveiwing of user'''
    gffs_f = []
    gffs_r = []
    fh = open(args_sc.main_gff)
    for entry in Gff3Parser().entries(fh):
        if entry.strand == "+":
            gffs_f.append(entry)
        else:
            gffs_r.append(entry)
    gffs_f = sorted(gffs_f, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    gffs_r = sorted(gffs_r, key=lambda k: (k.seq_id, k.start, k.end, k.strand))
    out_f = open(forward_file, "w")
    print_batch(args_sc, out_f, "+", libs["ft"],
                libs["fn"], libs["ff"], strain)
    out_r = open(reverse_file, "w")
    print_batch(args_sc, out_r, "-", libs["rt"],
                libs["rn"], libs["rf"], strain)
    seq = get_length(args_sc.fasta)
    gen_batch(libs["ft"], libs["fn"], libs["ff"], "+", gffs_f, out_f, seq)
    gen_batch(libs["rt"], libs["rn"], libs["rf"], "-", gffs_r, out_r, seq)
    fh.close()
    out_f.close()
    out_r.close()
