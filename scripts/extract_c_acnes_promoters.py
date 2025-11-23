#!/usr/bin/env python3
"""
extract_c_acnes_promoters.py

Input (expected files in current directory):
 - GCF_000008345.1_ASM834v1_genomic.fna    (genome FASTA)
 - GCF_000008345.1_ASM834v1_genomic.gff    (GFF3 file)

Output:
 - C_acnes_promoters_300bp.fasta

Defaults: UPSTREAM = 300 (bp). Change via --upstream argument.
"""

import sys
import argparse
from pathlib import Path
from Bio import SeqIO
from BCBio import GFF

def load_genome(fasta_path):
    seqdict = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        seqdict[rec.id] = rec.seq
    return seqdict

def parse_gff(gff_path):
    # we will parse features using BCBio.GFF to get feature locations precisely
    regions = []
    with open(gff_path) as g:
        for rec in GFF.parse(g):
            chrom = rec.id
            for feature in rec.features:
                # flatten nested features
                # capture CDS / gene entries
                if feature.type in ("CDS", "gene"):
                    attrs = feature.qualifiers
                    start = int(feature.location.start) + 1  # Biopython uses 0-based locations internally for SeqFeature; adjust if needed
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    # pick identifiers
                    locus = attrs.get("locus_tag", attrs.get("ID", [""]))[0]
                    name = attrs.get("Name", attrs.get("gene", [""]))[0] if attrs else ""
                    regions.append(dict(chrom=chrom, start=start, end=end, strand=strand, locus=locus, name=name))
    return regions

def build_occupied_intervals(regions):
    # intervals per chrom to check overlaps quickly
    from bisect import bisect_left, bisect_right
    iv = {}
    for r in regions:
        iv.setdefault(r['chrom'], []).append((r['start'], r['end']))
    # sort
    for chrom in iv:
        iv[chrom].sort()
    return iv

def overlaps_any(chrom, a, b, intervals):
    # naive check (regions list small for bacterial genome)
    for s,e in intervals.get(chrom, []):
        if not (b < s or a > e):
            return True
    return False

def extract_promoters(genome_fa, regions, upstream=300, output_fa="C_acnes_promoters_300bp.fasta"):
    from Bio.Seq import Seq
    with open(output_fa, "w") as out:
        seqdict = load_genome(genome_fa)
        intervals = build_occupied_intervals(regions)
        count = 0
        for r in regions:
            chrom = r['chrom']
            if chrom not in seqdict:
                # try trimming header differences: sometimes chrom id differs (use first seq)
                chrom = list(seqdict.keys())[0]
            seq = seqdict[chrom]
            if r['strand'] == 1 or r['strand'] == None:
                # positive strand -> upstream is before start
                a = max(1, r['start'] - upstream)
                b = r['start'] - 1
                # trim if overlaps any other feature (we avoid including cds of upstream gene)
                # naive: if overlap, reduce a to the end of last upstream feature +1
                # for simplicity here, just shorten to avoid negative region
            else:
                # negative strand -> upstream is after end (downstream on genome coords)
                a = r['end'] + 1
                b = min(len(seq), r['end'] + upstream)
            if a > b:
                # no space to extract
                continue
            subseq = seq[a-1:b]  # Python 0-based
            if r['strand'] == -1:
                subseq = subseq.reverse_complement()
            header = f">{r.get('locus','.') }|{r.get('name','.') }|{chrom}:{a}-{b}({ '+' if r['strand'] in (1,None) else '-' })"
            out.write(header + "\n")
            out.write(str(subseq) + "\n")
            count += 1
    print(f"Wrote {count} promoter sequences to {output_fa}")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--gff", default="GCF_000008345.1_ASM834v1_genomic.gff", help="GFF3 file")
    p.add_argument("--fasta", default="GCF_000008345.1_ASM834v1_genomic.fna", help="Genome FASTA")
    p.add_argument("--upstream", type=int, default=300, help="Upstream bp to extract")
    p.add_argument("--out", default="C_acnes_promoters_300bp.fasta", help="Output FASTA")
    args = p.parse_args()

    # dependency note: requires biopython and bcbio-gff
    # install via: pip install biopython bcbio-gff
    genome = Path(args.fasta)
    gff = Path(args.gff)
    if not genome.exists() or not gff.exists():
        print("Error: genome FASTA or GFF not found. Make sure filenames match.")
        sys.exit(1)

    regions = parse_gff(str(gff))
    extract_promoters(str(genome), regions, upstream=args.upstream, output_fa=args.out)
