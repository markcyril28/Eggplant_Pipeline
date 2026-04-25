#!/usr/bin/env python3
import csv

input_file = "plantCARE_output_PlantCARE_5913.tab"         # your original file
output_file = "plantCARE_output_PlantCARE_5913_reformatted.tab"

with open(input_file, "r", encoding="utf-8") as fin, open(output_file, "w", newline="", encoding="utf-8") as fout:
    reader = csv.reader(fin, delimiter="\t")
    writer = csv.writer(fout, delimiter="\t")
    for row in reader:
        if len(row) < 6:
            continue  # skip incomplete lines
        gene_id = row[0].strip()
        element = row[1].strip()
        try:
            start = int(row[3])
            length = int(row[4])
            end = start + length - 1
        except ValueError:
            continue
        strand = row[5].strip()
        writer.writerow([gene_id, element, start, end, strand])
