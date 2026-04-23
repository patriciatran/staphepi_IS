import os
from Bio import SeqIO
from pathlib import Path
import csv

# === Parameters ===
gbff_root = Path("ncbi_dataset/data")
output_dir = Path("dnaA_to_rlmH_fasta_per_assembly_1Mbp")
output_dir.mkdir(parents=True, exist_ok=True)

summary_path = output_dir / "dnaA_rlmH_coordinates_1Mbp.tsv"
downstream_buffer = 1_000_000
start_gene = "dnaA"
end_gene = "rlmH"

summary_rows = []

for gbff_file in gbff_root.glob("*/*.gbff"):
    with open(gbff_file, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            if "plasmid" in record.description.lower():
                continue

            dnaA_feature = None
            rlmH_feature = None

            for feature in record.features:
                if feature.type != "gene":
                    continue
                gene_names = feature.qualifiers.get("gene", [])
                if start_gene in gene_names:
                    dnaA_feature = feature
                elif end_gene in gene_names:
                    rlmH_feature = feature

            if not dnaA_feature or not rlmH_feature:
                print(f"Skipping {gbff_file.name}: missing dnaA or rlmH")
                continue

            dnaA_start = int(dnaA_feature.location.start)
            rlmH_strand = rlmH_feature.location.strand
            assembly = gbff_file.parts[-2]

            if rlmH_strand == 1:
                rlmH_start = int(rlmH_feature.location.start)
                rlmH_end = int(rlmH_feature.location.end)
                region_start = dnaA_start
                region_end = min(len(record.seq), rlmH_end + downstream_buffer)
                if region_end <= region_start:
                    continue
                subseq = record.seq[region_start:region_end]
                rlmH_rel_start = rlmH_start - dnaA_start
                rlmH_rel_end = rlmH_end - dnaA_start
            elif rlmH_strand == -1:
                rlmH_start = int(rlmH_feature.location.start)
                rlmH_end = int(rlmH_feature.location.end)
                region_start = max(0, rlmH_start - downstream_buffer)
                region_end = dnaA_start
                if region_end <= region_start:
                    continue
                subseq = record.seq[region_start:region_end].reverse_complement()
                # Reverse-complemented region: compute relative positions from end
                rlmH_rel_start = region_end - rlmH_end
                rlmH_rel_end = region_end - rlmH_start
            else:
                print(f"Skipping {gbff_file.name}: unknown rlmH strand")
                continue

            fasta_id = f"{assembly}_{record.id}_dnaA0_rlmHplus1Mbp"
            fasta_path = output_dir / f"{fasta_id}.fasta"
            with open(fasta_path, "w") as out_fasta:
                out_fasta.write(f">{fasta_id}\n")
                for i in range(0, len(subseq), 60):
                    out_fasta.write(str(subseq[i:i+60]) + "\n")

            summary_rows.append({
                "Assembly": assembly,
                "Contig_ID": record.id,
                "Strand": "+" if rlmH_strand == 1 else "-",
                "Extracted_Length": len(subseq),
                "rlmH_relative_start": rlmH_rel_start,
                "rlmH_relative_end": rlmH_rel_end
            })

# Write summary table
with open(summary_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary_rows[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"\nCompleted. Summary table: {summary_path}")
