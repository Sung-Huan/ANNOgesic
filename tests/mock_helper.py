from mock_gff3 import Create_generator

def convert_dict(line_list):
    datas = {}
    for data in line_list:
        datas[data] = data
    return datas

def gen_file(out_file, content):
    with open(out_file, "w") as fh:
        fh.write(content)

def import_data(filename):
    datas = []
    with open(filename) as fh:
        for line in fh:
            line = line.rstrip()
            datas.append(line)
    return datas

def extract_info(out_file, type_):
    datas = []
    attributes = []
    if type_ == "file":
        with open(out_file) as fh:
            for line in fh:
                line = line.rstrip()
                if len(line):
                    attributes.append(line.split("\t")[-1].split(";"))
                    datas.append("\t".join(line.split("\t")[0:-1]))
    else:
        for line in out_file.split("\n"):
            line = line.rstrip()
            if len(line):
                attributes.append(line.split("\t")[-1].split(";"))
                datas.append("\t".join(line.split("\t")[0:-1]))
    
    return datas, attributes

def read_dict(num, gff, attributes):
    gffs = []
    for index in range(0, num):
        gffs.append(Create_generator(gff[index], attributes[index], "gff"))
    return gffs
