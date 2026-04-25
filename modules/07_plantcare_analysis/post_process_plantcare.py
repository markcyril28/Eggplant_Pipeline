#!/usr/bin/env python3
import csv
import os
import argparse
import tarfile
import tempfile
import re


def derive_output_stem(input_file):
    """Build a stable output stem from .tab or .tar.gz input filenames."""
    base = os.path.basename(input_file)
    if base.endswith('.tar.gz'):
        return base[:-7]
    if base.endswith('.tab'):
        return base[:-4]
    return os.path.splitext(base)[0]


def derive_sequence_id(raw_sequence_id, output_stem):
    """Prefer gene ID from archive/filename (e.g., SMEL5_01g008730), fallback to raw PlantCARE ID."""
    stem = output_stem or ""

    # Most reliable match for this pipeline's gene naming convention.
    smel_match = re.search(r"(SMEL\d+_\d+g\d+)", stem)
    if smel_match:
        return smel_match.group(1)

    # Generic fallback: PlantCARE_<num>_<gene_token>_plantCARE
    generic_match = re.search(r"PlantCARE_\d+_([A-Za-z0-9_\-.]+?)_plantCARE(?:__|$)", stem)
    if generic_match:
        return generic_match.group(1)

    return raw_sequence_id


def process_plantcare_file(input_file, output_dir, output_stem=None):
    """
    Processes a single PlantCARE tab file and creates two post-processed versions.
    Uses gene ID parsed from filename when available; falls back to first-column PlantCARE ID.
    """
    if output_stem is None:
        output_stem = derive_output_stem(input_file)
    tbtools_output_file = os.path.join(output_dir, f"{output_stem}_tbtools.tab")
    heatmap_output_file = os.path.join(output_dir, f"{output_stem}_heatmap.tab")

    with open(input_file, "r", encoding="utf-8") as fin, \
         open(tbtools_output_file, "w", newline="", encoding="utf-8") as fout_tb, \
         open(heatmap_output_file, "w", newline="", encoding="utf-8") as fout_heat:

        reader = csv.reader(fin, delimiter="\t")
        writer_tb = csv.writer(fout_tb, delimiter="\t")
        writer_heat = csv.writer(fout_heat, delimiter="\t")

        # Write header for heatmap file
        writer_heat.writerow(['Sequence_ID', 'Motif_Name', 'Motif_Sequence', 'Position', 'Length', 'Strand', 'Organism', 'Function'])

        for row in reader:
            if len(row) < 6:
                continue  # skip incomplete lines

            raw_sequence_id = row[0].strip()
            gene_id = derive_sequence_id(raw_sequence_id, output_stem)
            motif_name = row[1].strip()
            if not motif_name:
                continue # skip rows with no motif name

            # TBTools format
            try:
                start = int(row[3])
                length = int(row[4])
                end = start + length - 1
            except (ValueError, IndexError):
                continue
            strand = row[5].strip()
            writer_tb.writerow([gene_id, motif_name, start, end, strand])

            # Heatmap format - ensure we have all 8 columns
            heatmap_row = [gene_id] + row[1:8] if len(row) >= 8 else [gene_id] + row[1:] + [''] * (8 - len(row))
            writer_heat.writerow(heatmap_row)

def main():
    """Main function to run the post-processing."""
    parser = argparse.ArgumentParser(description='Post-process PlantCARE tab files.')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory with raw PlantCARE tab files.')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for post-processed files.')
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' not found.")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    # Recursively find .tab files and PlantCARE .tar.gz archives (skip hidden dirs)
    processed = 0
    for dirpath, dirnames, filenames in os.walk(args.input_dir):
        dirnames[:] = [d for d in dirnames if not d.startswith('_')]
        for filename in filenames:
            if filename.endswith(".tab"):
                input_file = os.path.join(dirpath, filename)
                process_plantcare_file(input_file, args.output_dir)
                print(f"Processed '{filename}' from {os.path.relpath(dirpath, args.input_dir)}")
                processed += 1
            elif filename.endswith(".tar.gz"):
                archive_path = os.path.join(dirpath, filename)
                archive_stem = derive_output_stem(archive_path)
                try:
                    with tarfile.open(archive_path, "r:gz") as tar:
                        tab_members = [m for m in tar.getmembers() if m.isfile() and m.name.endswith(".tab")]
                        if not tab_members:
                            print(f"Warning: No .tab file found in archive '{filename}'")
                            continue

                        with tempfile.TemporaryDirectory() as temp_dir:
                            for member in tab_members:
                                extracted_name = os.path.basename(member.name)
                                extracted_path = os.path.join(temp_dir, extracted_name)
                                extracted_file = tar.extractfile(member)
                                if extracted_file is None:
                                    continue
                                with open(extracted_path, "wb") as fout:
                                    fout.write(extracted_file.read())

                                member_stem = derive_output_stem(extracted_name)
                                out_stem = f"{archive_stem}__{member_stem}"
                                process_plantcare_file(extracted_path, args.output_dir, output_stem=out_stem)
                                processed += 1

                        print(f"Processed archive '{filename}' from {os.path.relpath(dirpath, args.input_dir)}")
                except tarfile.TarError:
                    print(f"Warning: Could not read archive '{filename}'")

    if processed == 0:
        print(f"Warning: No .tab files or .tar.gz archives found in '{args.input_dir}' (searched recursively).")

if __name__ == '__main__':
    main()
