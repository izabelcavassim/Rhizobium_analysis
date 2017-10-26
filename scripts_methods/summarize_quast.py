import glob

base_dir = "/home/agb/data/rhizo/full_dataset/quast"

methods = ["SPAdes", "Abyss", "a5", "CISA"]
ignore = ["3211-8", "3260-63", "3261-65", "3262-66", "3271-76", "3278-83", "3281-86", "3282-87", "3290-95", "3293-98",
          "3339-136B", "3381-156"]

contigs_0bp = {m: [] for m in methods}
contigs_500bp = {m: [] for m in methods}
total_length_0bp = {m: [] for m in methods}
total_length_500bp = {m: [] for m in methods}
largest_contig = {m: [] for m in methods}
gc_content = {m: [] for m in methods}
n50 = {m: [] for m in methods}
n75 = {m: [] for m in methods}
l50 = {m: [] for m in methods}
l75 = {m: [] for m in methods}

for tsv_path in glob.glob("{}/*.tsv".format(base_dir)):
    if tsv_path[(tsv_path.rindex("/") + 1):-4] in ignore:
        continue

    with open(tsv_path) as file_h:
        lines = [line.split("\t") for line in file_h.readlines()]
        for idx, method in enumerate(methods):
            contigs_0bp[method].append(int(lines[1][idx + 1]))
            contigs_500bp[method].append(int(lines[5][idx + 1]))
            total_length_0bp[method].append(int(lines[3][idx + 1]))
            total_length_500bp[method].append(int(lines[7][idx + 1]))
            largest_contig[method].append(int(lines[6][idx + 1]))
            gc_content[method].append(float(lines[8][idx + 1]))
            n50[method].append(int(lines[9][idx + 1]))
            n75[method].append(int(lines[10][idx + 1]))
            l50[method].append(int(lines[11][idx + 1]))
            l75[method].append(int(lines[12][idx + 1]))


def avg(vals): return float(sum(vals)) / len(vals)


def med(vals): return sorted(vals)[int(len(vals) / 2)]


def print_line(label, vals):
    print("{}\t{}".format(label, "\t".join(str(round(avg(vals[method]), 2)) for method in methods)))


print("Method\t{}".format("\t".join(method for method in methods)))
print_line("# contigs (>= 0 bp)", contigs_0bp)
print_line("# contigs (>= 500 bp)", contigs_500bp)
print_line("Total length (>= 0 bp)", total_length_0bp)
print_line("Total length (>= 500 bp)", total_length_500bp)
print_line("largest contig", largest_contig)
print_line("GC (%)", gc_content)
print_line("N50", n50)
print_line("N75", n75)
print_line("L50", l50)
print_line("L75", l75)
