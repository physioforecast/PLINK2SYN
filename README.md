# PLINK2SYN

Synthetic data generation and GWAS validation pipeline using PLINK2 and 1000 Genomes Project data.

## Overview

PLINK2SYN is a genomics research pipeline for generating synthetic genotype-phenotype datasets with known genetic effects. This allows validation and testing of GWAS (Genome-Wide Association Study) analysis methods by comparing results against ground truth.

The pipeline:
1. Processes real genetic data from the 1000 Genomes Project Phase 3
2. Generates synthetic samples by adding controlled noise to genotypes
3. Creates phenotypes with known heritability and genetic architecture
4. Runs GWAS using both PLINK2 and R implementations
5. Enables comparison of results to validate statistical methods

## Features

- **Data Augmentation**: Generate ~32K synthetic samples from real genotypes with controlled Gaussian noise
- **Flexible Missingness**: Test scenarios with and without missing genotype data
- **Known Ground Truth**: Save true causal variants for validation
- **Dual Implementation**: Run GWAS in both PLINK2 and R for cross-validation
- **Memory Efficient**: Chunk-based processing for large datasets
- **Configurable Parameters**: Customizable heritability, noise levels, missing data rates

## Requirements

### Software Dependencies
- **PLINK2** (v2.0 or later) - Download from [PLINK2 website](https://www.cog-genomics.org/plink/2.0/)
- **Python 3.8+** with packages:
  - numpy
  - pandas
  - scipy
  - tqdm
- **R 4.0+** with packages:
  - data.table
  - logistf
  - jsonlite
- **Bash** (for orchestration scripts)

### Data Requirements
- 1000 Genomes Project Phase 3 data in PLINK2 format (pgen/pvar/psam)
  - Download from [PLINK2 Resources - 1000 Genomes Phase 3](https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg)
  - Files should be named: `all_hg38.pgen`, `all_hg38.pvar.zst`, `all_hg38.psam`

## Setup

### 1. Clone Repository
```bash
git clone <repository-url>
cd PLINK2SYN
```

### 2. Set Up Python Environment
```bash
python3 -m venv PLINK2SYN-venv
source PLINK2SYN-venv/bin/activate
pip install -r requirements.txt
```

### 3. Install PLINK2
Download the appropriate PLINK2 binary for your platform and ensure it's in your PATH:
```bash
# Example for macOS ARM64
wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_mac_arm64_<version>.zip
unzip plink2_mac_arm64_<version>.zip
# Add to PATH or create symlink
```

### 4. Download Source Data
Place 1000 Genomes Project Phase 3 data in the `sourcedata/` directory:
```
sourcedata/
├── all_hg38.pgen
├── all_hg38.pvar.zst
└── all_hg38.psam
```

### 5. Configure Paths
Create a `config.sh` file with your derivative data path:
```bash
export derivativepath="/path/to/your/output/directory/"
```

## Usage

### Generate Base Synthetic Data

The master script processes both scenarios (with and without missing data):

```bash
./master_makebase.sh
```

This will:
1. Filter to 50K SNPs using PLINK2
2. Export to additive variance format
3. Generate ~32K synthetic samples in chunks
4. Combine chunks into final datasets

### Generate Phenotypes and Run GWAS

```bash
./master_makeglm.sh
```

This will:
1. Generate synthetic phenotypes with known genetic effects
2. Create both quantitative and binary phenotypes
3. Save ground truth causal variants
4. Produce phenotype/covariate files for GWAS

### Run GWAS Analysis

**Using PLINK2:**
```bash
source code/PLINK2R.sh
```

**Using R:**
```bash
Rscript code/RUNGWAS_R.R
```

## Workflow Details

### Step 1: Data Filtering and Export
- Selects biallelic SNPs only
- Thins to 50,000 random variants
- Optionally filters by genotyping rate (--geno 0 for nomiss scenario)
- Exports to PLINK2 Av (additive variance) format

### Step 2: Data Augmentation
**Script:** `code/augment_1kGP3.py`

Parameters:
- `--nsamples`: Number of synthetic samples (default: 32017)
- `--mean`: Mean of Gaussian noise (default: 0)
- `--std_dev`: Standard deviation of noise (default: 0.1)
- `--miss`: Enable additional missingness
- `--p`: Probability of missing values (default: 0.1)

Example:
```bash
python code/augment_1kGP3.py \
  --fname derivatives/data.traw \
  --fname2 sourcedata/all_hg38.psam \
  --nsamples 8000 \
  --miss --p 0.1
```

### Step 3: Phenotype Generation
**Script:** `code/gen_gwas.py`

Parameters:
- `--n_vars`: Number of covariates (default: 5)
- `--h2`: Heritability (default: 0.5)
- `--theta`: Case prevalence threshold (default: 0.1)

Example:
```bash
python code/gen_gwas.py derivatives/data.traw --h2 0.5 --theta 0.1
```

Outputs:
- `*_phenocov.csv`: Phenotypes and covariates
- `*_known_g_effects.csv`: Ground truth causal variants

### Step 4: GWAS Analysis
Run association testing using:
- PLINK2's `--glm` command
- R's logistic/linear regression (including Firth regression)

Compare results to validate methods against known causal variants.

## File Structure

```
PLINK2SYN/
├── README.md                      # This file
├── requirements.txt               # Python dependencies
├── .gitignore                     # Git exclusions
├── config.sh                      # User-specific paths (not tracked)
├── master_makebase.sh             # Main pipeline: base data generation
├── master_makeglm.sh              # Main pipeline: phenotype generation
├── code/
│   ├── augment_1kGP3.py          # Synthetic sample generation
│   ├── gen_gwas.py               # Phenotype generation
│   ├── syntest_toAv.sh           # PLINK2 processing orchestration
│   ├── PLINK2R.sh                # PLINK2 GWAS wrapper
│   ├── PLINK2R_subset.sh         # Subset GWAS analysis
│   └── RUNGWAS_R.R               # R-based GWAS implementation
├── sourcedata/                    # 1000 Genomes data (not tracked)
│   ├── all_hg38.pgen
│   ├── all_hg38.pvar.zst
│   └── all_hg38.psam
├── test_data/                     # Generated datasets (not tracked)
└── PLINK2SYN-venv/                # Python environment (not tracked)
```

## Testing Scenarios

### 1. No Missing Data (nomiss)
- Filters out any variants with missing genotypes
- Uses `--geno 0` in PLINK2
- Ideal for testing standard logistic/linear regression

### 2. With Missing Data (yesmiss)
- Preserves original missingness
- Adds ~10% additional random missingness
- Tests handling of missing data in statistical models

## Output Files

### Genotype Files
- `*.pgen`, `*.pvar`, `*.psam` - PLINK2 binary format
- `*.traw` - Transposed raw format (for Python processing)

### Phenotype Files
- `*_phenocov.csv` - Contains:
  - `FID`, `IID` - Family and individual IDs
  - `y` - Quantitative phenotype
  - `ybool` - Binary phenotype (1=control, 2=case)
  - `COV_1` through `COV_N` - Correlated covariates

### Ground Truth Files
- `*_known_g_effects.csv` - Contains:
  - `variant IDs` - SNP identifiers
  - `x_values` - Binary indicator (1=causal, 0=null)

## Parameters and Tuning

### Heritability Model
The phenotype generation uses:
- σ²_genetic = h² (default: 0.5)
- σ²_covariate = 0.1
- σ²_environmental = 1 - h² - 0.1

### Noise Model
Genotype augmentation:
- Adds N(0, 0.1) noise to dosages
- Clips values to [0, 2] range
- Preserves allele dosage distribution

### Causal Variants
- Selected with probability 0.01 (~500 causal out of 50K)
- Effect sizes standardized through z-scoring
- Contributions weighted by heritability parameter

## What's Not in Git

This repository uses a code-only approach. The following are excluded via `.gitignore`:

- **Data files** (`/sourcedata/`, `/test_data/`) - Multi-GB genetic datasets
- **PLINK2 binaries** - Platform-specific executables
- **Python environment** (`PLINK2SYN-venv/`) - Virtual environment and packages
- **Configuration** (`config.sh`) - Machine-specific paths
- **Results** (`PLINK2SYN_*`) - Generated analysis outputs

Users must:
1. Download 1000 Genomes data separately
2. Install PLINK2 for their platform
3. Create their own Python environment
4. Configure paths for their system
5. Generate datasets by running the pipeline

## Validation Workflow

1. Generate synthetic data with known causal variants
2. Run GWAS using PLINK2 and/or R
3. Load results and ground truth:
   ```python
   import pandas as pd
   gwas_results = pd.read_csv('gwas_output.glm.logistic', sep='\t')
   true_effects = pd.read_csv('*_known_g_effects.csv')
   ```
4. Compare discovered associations to known causal variants
5. Calculate metrics (sensitivity, specificity, PPV, etc.)

## Memory Considerations

The pipeline uses chunk-based processing to handle large datasets:
- Default chunk size: 8,000 samples
- Processes ~32K samples in 4-5 chunks
- Intermediate files are combined after processing
- Monitor disk space for derivative files (~10-20GB)

## Troubleshooting

**"Neither python nor python3 found"**
- Activate the virtual environment: `source PLINK2SYN-venv/bin/activate`

**"PLINK2 not found"**
- Ensure PLINK2 is in your PATH or update scripts with full path

**Memory errors**
- Reduce `chunk_size` in `syntest_toAv.sh`
- Reduce `nsamples` parameter
- Ensure sufficient RAM (recommend 16GB+)

**Missing data warnings**
- Check that source data files are properly formatted
- Verify file paths in `config.sh`

<!-- ## Citation

If you use this pipeline in your research, please cite:
- PLINK2: Chang CC, Chow CC, Tellier LC, et al. (2015)
- 1000 Genomes Project: The 1000 Genomes Project Consortium (2015)

## License

[Specify your license here]

## Contact

[Add contact information or remove this section] -->
