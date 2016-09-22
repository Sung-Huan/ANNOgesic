def change_format(input_file, output_file):
    '''change the format of sRNA database'''
    num = 1
    out = open(output_file, "w")
    with open(input_file) as f_h:
        for line in f_h:
            line = line.strip()
            if line.startswith(">"):
                datas = line.split("|")
                if datas[0][1:] == "NA":
                    datas[0] = ">srn_" + str(num)
                    num += 1
                out.write("|".join(datas[:3]) + "\n")
            else:
                out.write(line + "\n")
    out.close()
