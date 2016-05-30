import argparse

__author__ = "Sung-Huan Yu <sung-huan.yu@uni-wuerzburg.de>"
__email__ = "sung-huan.yu@uni-wuerzburg.de"

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file",help="input file")
parser.add_argument("-o","--output_file",help="output file")
parser.add_argument("-r","--ratt_path",help="output file")
parser.add_argument("-p","--pagit_path",help="output file")
args = parser.parse_args()

def main():
    out = open(args.output_file, "w")
    with open(args.input_file) as fg:
        for line in fg:
            datas = line.split("=")
            if datas[0] == "RATT_HOME":
                out.write("RATT_HOME=" + args.ratt_path + "; export RATT_HOME\n")
            elif datas[0] == "NUCMER_PATH":
                out.write("NUCMER_PATH=" + args.pagit_path + "/bin/;\n")
            else:
                out.write(line)
if __name__ == "__main__":
    main()
