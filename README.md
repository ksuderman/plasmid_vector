# C. acnes Promoter Strength Prediction Project Guide

## Project Overview

**Goal**: Predict which promoter sequences will be strongest for plasmid vectors in C. acnes chassis, with focus on improving weak P3 promoter from S. aureus through -10/-35 region modifications.

**Background**: Previous research showed that altering -10 and -35 promoter regions can increase GFP activity in plasmids. This project aims to use machine learning to predict optimal promoter sequences for C. acnes.

## Key Research Findings

### Promoter-GPT Analysis
**Important**: Promoter-GPT is NOT suitable for this project because:
- Designed for generating new sequences, not predicting strength of existing ones
- 0.43M parameter transformer that creates novel 200bp promoter sequences
- Cannot analyze or rank existing promoters for strength
- Uses k-mer tokenization (k=3) for sequence generation

### Better ML Approaches for Promoter Strength Prediction

#### 1. msBERT-Promoter (2024) - Most Promising
- Two-stage predictor: identifies promoters AND predicts their strengths
- Uses BERT-based architecture fine-tuned for DNA sequences
- Incorporates multi-scale sequence information
- Published in BMC Biology, 2024

#### 2. XGBoost-Based Models (2024)
- Excellent performance metrics: R² = 0.88, correlation = 0.94
- Specifically designed for predicting strength of artificially designed promoters
- Could be adapted for C. acnes sequences
- Published in ACS Synthetic Biology

#### 3. MLDSPP (2024)
- Machine Learning and Duplex Stability based Promoter prediction
- Uses DNA structural properties with explainable AI
- Tested on 12 diverse bacterial genomes
- Published in Journal of Chemical Information and Modeling

### Research Gap Identified
**No existing models are specifically trained on C. acnes promoters** - this creates an opportunity for novel research in this important skin microbiome organism.

## Technical Background

### Bacterial Promoter Elements
- **-35 region**: Consensus sequence TTGACA
- **-10 region**: Consensus sequence TATAAT (Pribnow box)
- **Spacer region**: Typically 16-18 bp between -35 and -10
- **UP element**: Additional regulatory sequence upstream of -35

### -10/-35 Region Modification Strategy
Based on classical bacterial promoter research:
- Closer matches to consensus sequences generally increase strength
- -10 element (T-12A-11T-10A-9A-8T-7) has extensive protein-DNA interactions
- Base-specific interactions primarily with A-11 and T-7 positions

## Recommended Implementation Strategy

### Phase 1: Data Collection & Preparation
1. Compile C. acnes promoter sequences (target: 50-100 sequences)
2. Include weak P3 promoter from S. aureus as baseline
3. Format as FASTA sequences (50-200bp upstream of start codons)
4. Gather any available experimental strength data (GFP, qPCR, etc.)

### Phase 2: Model Selection & Implementation
**Option A**: Transfer Learning with msBERT-Promoter
- Use pre-trained BERT model fine-tuned for promoters
- Adapt for C. acnes-specific sequences
- Better for complex sequence patterns

**Option B**: Feature-Based XGBoost Model
- Extract sequence features (GC content, motif presence, etc.)
- Train XGBoost regressor on promoter strength
- More interpretable results

### Phase 3: -10/-35 Region Optimization
1. Use trained model to score current P3 promoter
2. Generate systematic mutations in -10/-35 regions
3. Score all variants with model
4. Select top candidates for experimental validation

## Computational Requirements

### Software Environment
- Python 3.8+
- PyTorch or TensorFlow for deep learning models
- scikit-learn for traditional ML models
- BioPython for sequence manipulation
- pandas/numpy for data handling

### Hardware Recommendations
- Minimum: 8GB RAM, modern CPU
- Optimal: GPU support for transformer models
- Storage: ~1GB for model weights and data

## Next Steps Priority List

1. **Immediate**: Set up computational environment
2. **Data Collection**: Gather C. acnes promoter sequences
3. **Model Implementation**: Choose and implement prediction approach
4. **Validation**: Test against known strong/weak promoters
5. **Design**: Create -10/-35 region modifications
6. **Experimental**: Validate top predictions in lab

## Research Sources

### Key Papers (2024)
- [MLDSPP: Bacterial Promoter Prediction Tool](https://pubs.acs.org/doi/10.1021/acs.jcim.3c02017)
- [Precise Prediction of Promoter Strength](https://pubs.acs.org/doi/abs/10.1021/acssynbio.1c00117)
- [msBERT-Promoter](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-01923-z)
- [Combining diffusion and transformer models](https://journals.asm.org/doi/10.1128/msystems.00183-25)
- [Engineered probiotic Lactobacillus plantarum](https://journals.asm.org/doi/10.1128/spectrum.01829-23)

### Additional Resources
- [Promoter-GPT Blog](https://huggingface.co/blog/hugging-science/promoter-gpt)
- [Online Analysis Tools - Promoters](https://molbiol-tools.ca/Promoters.htm)
- [PromoterLCNN](https://pmc.ncbi.nlm.nih.gov/articles/PMC9325283/)

## Questions for Project Planning

1. How many C. acnes promoter sequences are available?
2. Is there existing experimental data on promoter strengths?
3. What is the preferred programming experience level?
4. What computational resources are available?
5. Timeline for experimental validation?

---

*Document created: November 2024*
*Status: Planning phase - ready for computational setup*