import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot
import pylab
import sys

def Plot_mountain_plot(input_file, output_name)
    poss = []
    values = []
    check = 0
    pre_check = 0
    fh = open(args.input_file, "r");
    while True:
        line = fh.readline()
        line = line.rstrip()
        if not line:
            matplotlib.pyplot.figure(1)
            matplotlib.pyplot.subplot(212)
            matplotlib.pyplot.xlabel('Nucleotide position')
            matplotlib.pyplot.ylabel('Entropy')
            matplotlib.pyplot.plot(values, color="black")
            matplotlib.pyplot.savefig(args.output_name, format='pdf')
            break
        elif line == "&":
            line = fh.readline()
            line = line.rstrip()
            check += 1
        else:
            poss.append(float(line[0:4].replace(" ", "")))
            values.append(float(line[5:].replace(" ","")))
        if check != pre_check:
            pre_check = check
            if check == 1:
                matplotlib.pyplot.figure(1)
                matplotlib.pyplot.subplot(211)
                matplotlib.pyplot.ylabel('Number of enclosing nucleotides\nor\nMin free energy structure', fontsize=10, multialignment='left')
                matplotlib.pyplot.plot(values, label='pair probabilities')
                values = []
                poss = []
            elif check == 2:
                matplotlib.pyplot.plot(values, label='mfe structure')
                matplotlib.pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.)
                values = []
                poss = []
    fh.close()
