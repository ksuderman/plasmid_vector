#!/usr/bin/env bash
set -eu

# Example: download GenBank and GFF3 + chromosome FASTA from NCBI Assembly FTP (adjust if needed)
# Replace GCF_000008345.1_ASM834v1 with another assembly path if desired.

BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/345/GCF_000008345.1_ASM834v1"
wget "$BASE_URL/GCF_000008345.1_ASM834v1_genomic.fna.gz"
wget "$BASE_URL/GCF_000008345.1_ASM834v1_genomic.gff.gz"
wget "$BASE_URL/GCF_000008345.1_ASM834v1_genomic.gbff.gz"  # optional GenBank flat file

gunzip *.gz
