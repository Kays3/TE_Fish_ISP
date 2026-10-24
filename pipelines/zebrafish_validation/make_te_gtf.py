"""UCSC rmsk.txt.gz -> TEtranscripts-format TE GTF with Ensembl chromosome names.
gene_id = repName (subfamily), transcript_id = repName_dupN, family_id = repFamily,
class_id = repClass. Only chromosomes present in the Ensembl genome are kept."""
import gzip, sys, collections
rmsk, fasta_index, out = sys.argv[1:4]
ens = {l.split("\t")[0] for l in open(fasta_index)}
dup = collections.Counter(); kept = skipped = 0
with gzip.open(rmsk, "rt") as fh, open(out, "w") as o:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        chrom = f[5][3:] if f[5].startswith("chr") else f[5]
        chrom = "MT" if chrom == "M" else chrom
        if chrom.startswith("Un_"):                 # chrUn_KN149680v1 -> KN149680.1
            acc = chrom[3:]
            if "v" in acc:
                base, ver = acc.rsplit("v", 1)
                chrom = f"{base}.{ver}"
        if chrom not in ens:
            skipped += 1; continue
        name, cls, fam = f[10], f[11], f[12]
        dup[name] += 1
        o.write(f'{chrom}\trmsk\texon\t{int(f[6]) + 1}\t{f[7]}\t{f[1]}\t{f[9]}\t.\t'
                f'gene_id "{name}"; transcript_id "{name}_dup{dup[name]}"; family_id "{fam}"; class_id "{cls}";\n')
        kept += 1
print(f"kept {kept} loci ({len(dup)} subfamilies); skipped {skipped} on non-Ensembl contigs")
