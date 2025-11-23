# Data Collection Guide for C. acnes Promoter Analysis

## What You Need to Collect

### 1. P3 Promoter from S. aureus
**Source**: Look up the paper you referenced (https://doi.org/10.1128/spectrum.01829-23)
- Find the P3 promoter sequence used in that study
- Should be the weak promoter they improved by modifying -10/-35 regions
- If not available in paper, check supplementary materials

### 2. Strong C. acnes Promoter Sequences

#### Option A: Literature Search
Search for papers on C. acnes gene expression, promoter analysis, or recombinant systems:
- PubMed search: "Cutibacterium acnes promoter" or "C. acnes gene expression"
- Look for papers with promoter sequences or expression data
- Check supplementary materials for sequence data

#### Option B: Database Search
**NCBI RefSeq Genome Database**:
1. Go to https://www.ncbi.nlm.nih.gov/genome/
2. Search "Cutibacterium acnes"
3. Select a complete genome
4. Look for highly expressed genes (ribosomal proteins, housekeeping genes)
5. Extract ~100-200bp upstream of start codons

**Recommended C. acnes genes with likely strong promoters**:
- `rplL` (ribosomal protein L7/L12)
- `rpsA` (30S ribosomal protein S1)
- `groEL` (molecular chaperone)
- `recA` (DNA repair protein)
- `gyrA` (DNA gyrase)
- `dnaA` (chromosomal replication initiator)

#### Option C: PromoterBase Databases
- **RegulonDB**: http://regulondb.ccg.unam.mx/ (mostly E. coli but has bacterial promoter patterns)
- **DBTBS**: http://dbtbs.hgc.jp/ (Bacillus subtilis, similar Gram-positive)

### 3. Experimental Data (Optional but Helpful)
If you find papers with expression data:
- GFP reporter assay results
- qPCR expression levels
- Western blot quantification
- Any relative strength measurements

## How to Format Your Data

### FASTA Format Requirements:
```
>Promoter_Name_Description
SEQUENCE_HERE
```

### Example:
```
>rplL_promoter_strong
TTGACAGATCGATCGATCGATCGATCGATATAATAGCGATCGATCGATCGATC

>groEL_promoter_medium
TTGAAAGATCGATCGATCGATCGATCGATATAATAGCGATCGATCGATCGATC
```

### Important Notes:
1. **Sequence Length**: 50-200 bp upstream of start codon
2. **Include Context**: Make sure -35 and -10 regions are included
3. **Orientation**: 5' to 3' direction of the promoter
4. **Quality**: Use sequences from peer-reviewed sources

## Quick Start Steps

1. **Find P3 sequence**: Look up the spectrum.01829-23 paper
2. **Pick 5-10 C. acnes promoters**: Start with housekeeping genes
3. **Get sequences**: Extract from NCBI or literature
4. **Format as FASTA**: Use templates provided
5. **Test analysis**: Run sequence_analyzer.py

## Need Help?

If you're having trouble finding sequences:
1. Start with just P3 and 2-3 C. acnes promoters
2. Use example sequences to test the analysis pipeline
3. Gradually add more promoter sequences as you find them

## Resources

- **NCBI Genome Browser**: https://www.ncbi.nlm.nih.gov/genome/
- **PubMed**: https://pubmed.ncbi.nlm.nih.gov/
- **Paper you mentioned**: https://doi.org/10.1128/spectrum.01829-23