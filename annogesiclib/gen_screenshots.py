import os
import sys
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

def print_batch(out, strand, lib_t, lib_n, lib_f, fasta, main_gff, presentation,
                side_gffs, height, screenshot_folder, strain):
    out.write("new\n")
    out.write("genome {0}\n".format(os.path.join(os.getcwd(), fasta)))
    out.write("load {0}\n".format(os.path.join(os.getcwd(), main_gff)))
    gff = main_gff.split("/")
    out.write("{0} {1}\n".format(presentation, gff[-1]))
    for filename in side_gffs:
        out.write("load {0}\n".format(os.path.join(os.getcwd(), filename)))
        gff = filename.split("/")
        out.write("{0} {1}\n".format(presentation, gff[-1]))
    load_wigs(out, lib_t, lib_n, lib_f)
    out.write("maxPanelHeight {0}\n".format(height))
    if strand == "+":
        out.write("snapshotDirectory {0}\n".format(
                  os.path.join(os.getcwd(), screenshot_folder,
                  strain, "forward")))
    else:
        out.write("snapshotDirectory {0}\n".format(
                  os.path.join(os.getcwd(), screenshot_folder,
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

def gen_batch(lib_t, lib_n, lib_f, strand, gffs, out):
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
        out.write("goto {0}:{1}-{2}\n".format(
                  gff.seq_id, gff.start - 200, gff.end + 200))
        set_data_range(out, gff, wigs, strand)
        out.write("snapshot {0}:{1}-{2}.png\n".format(
                  gff.seq_id, gff.start, gff.end))

def gen_screenshot(main_gff, forward_file, reverse_file, screenshot_folder,
                   height, lib_ft, lib_fn, lib_rt, lib_rn, lib_ff, lib_rf,
                   fasta, side_gffs, presentation, strain):
    gffs_f = []
    gffs_r = []
    fh = open(main_gff)
    for entry in Gff3Parser().entries(fh):
        if entry.strand == "+":
            gffs_f.append(entry)
        else:
            gffs_r.append(entry)
    gffs_f = sorted(gffs_f, key=lambda k: (k.seq_id, k.start))
    gffs_r = sorted(gffs_r, key=lambda k: (k.seq_id, k.start))
    out_f = open(forward_file, "w")
    print_batch(out_f, "+", lib_ft, lib_fn, lib_ff, fasta, main_gff,
                presentation, side_gffs, height, screenshot_folder, strain)
    out_r = open(reverse_file, "w")
    print_batch(out_r, "-", lib_rt, lib_rn, lib_rf, fasta, main_gff,
                presentation, side_gffs, height, screenshot_folder, strain)
    gen_batch(lib_ft, lib_fn, lib_ff, "+", gffs_f, out_f)
    gen_batch(lib_rt, lib_rn, lib_rf, "-", gffs_r, out_r)
    fh.close()
    out_f.close()
    out_r.close()
