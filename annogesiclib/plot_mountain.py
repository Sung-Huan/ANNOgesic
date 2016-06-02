import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot
matplotlib.pyplot.style.use('ggplot')


def plot_mountain_plot(input_file, output_name):
    poss = []
    values = []
    check = 0
    pre_check = 0
    f_h = open(input_file, "r")
    while True:
        line = f_h.readline()
        line = line.rstrip()
        if not line:
            matplotlib.pyplot.figure(1)
            matplotlib.pyplot.subplot(212)
            matplotlib.pyplot.xlabel('Nucleotide position')
            matplotlib.pyplot.ylabel('Entropy')
            matplotlib.pyplot.plot(values, color="black")
            matplotlib.pyplot.savefig(output_name, format='pdf')
            break
        elif line == "&":
            line = f_h.readline()
            line = line.rstrip()
            check += 1
        else:
            poss.append(float(line[0:4].replace(" ", "")))
            values.append(float(line[5:].replace(" ", "")))
        if check != pre_check:
            pre_check = check
            if check == 1:
                matplotlib.pyplot.figure(1)
                matplotlib.pyplot.subplot(211)
                ylabel = ("Number of enclosing nucleotides\nor\n"
                          "Min free energy structure")
                matplotlib.pyplot.ylabel(
                    ylabel, fontsize=10, multialignment='left')
                matplotlib.pyplot.plot(values, label='pair probabilities')
                values = []
                poss = []
            elif check == 2:
                matplotlib.pyplot.plot(values, label='mfe structure')
                matplotlib.pyplot.legend(
                    bbox_to_anchor=(0., 1.02, 1., .102),
                    loc=3, ncol=2, mode="expand", borderaxespad=0.)
                values = []
                poss = []
    f_h.close()
    matplotlib.pyplot.cla()
    matplotlib.pyplot.clf()
