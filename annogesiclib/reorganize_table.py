import csv
import shutil
from annogesiclib.lib_reader import read_libs


def import_covers(row):
    cover_names = []
    covers = []
    for data in row.split("("):
        if ")" not in data:
            cover_names.append(data)
        else:
            covers.append(data.split(")")[0])
            if len(data.split(");")) == 2:
                cover_names.append(data.split(");")[-1])
    return cover_names, covers

def get_lib_name(libs):
    tracks = []
    double_tracks = []
    track_list = []
    for lib1 in libs:
        for lib2 in libs:
            if (lib1["cond"] == lib2["cond"]) and (
                    lib1["type"] == lib2["type"]) and (
                    lib1["rep"] == lib2["rep"]) and (
                    lib1["strand"] != lib2["strand"]):
                track = "/".join([lib1["name"], lib2["name"]])
                if track not in double_tracks:
                    double_tracks.append(track)
                    double_tracks.append("/".join([lib2["name"],
                                         lib1["name"]]))
                    tracks.append(track)
                    track_list.append([lib1["name"], lib2["name"]])
    return tracks, track_list

def reorganize_table(input_libs, wigs, cover_header, table_file):
    libs, texs = read_libs(input_libs, wigs)
    fh = open(table_file, "r")
    first = True
    headers = []
    tracks, track_list = get_lib_name(libs)
    out = open(table_file + "tmp", "w")
    for row in csv.reader(fh, delimiter='\t'):
        if first:
            detect = False
            header_num = 0
            for header in row:
                if header == cover_header:
                    index = header_num
                    detect = True
                header_num += 1
                if not detect:
                   headers.append(header)
                else:
                   detect = False
            first = False
            for track in tracks:
                headers.append("Avg_coverage:" + track)
            out.write("\t".join(headers) + "\n")
        else:
            if len(row) < (index + 1):
                cover_names = []
                covers = []
            else:
                cover_names, covers = import_covers(row[index])
            if len(row) == index + 1:
                row = row[:index]
            else:
                row = row[:index] + row[index + 1:]
            detects = ["Not_detect"] * len(tracks)
            for name, cover in zip(cover_names, covers):
                num_track = 0
                for track in track_list:
                    if name in track:
                        detects[num_track] = cover
                    num_track += 1
            out.write("\t".join(row + detects) + "\n")
    out.close()
    shutil.move(table_file + "tmp", table_file)
