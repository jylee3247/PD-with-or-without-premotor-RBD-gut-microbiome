# Introduction
This repository includes custom codes used for the analyses of this manuscripts.<br/>
The codes have been tested on Ubuntu20.04
- Title: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders: exploring the alpha-synuclein pathways
- Authors: Jae-Yun Lee, Sungyang Jo, Jihyun Lee, Sangjin Lee, Hyun Sik Kim, Jin-Woo Bae, Sun Ju Chung.
- Year: 2024

# Contents
* **16S_amplicon**: A custom script for 16S rRNA gene amplicon sequencing data analysis.
  - `Primer_trimming.py`: Used for PCR primer sequences from sequencing raw data. <br>
  - After primer trimming, the data were imported to [Qiime2 platform](https://qiime2.org/) for further analyses: denoising, taxonomy classification, filtering and phylogenetic tree construction.
 
* **Shotgun_metagenome**: Custom scripts for shotgun metagenome sequencing data analysis.
  - biobakery: Script for run [bioBakery](https://huttenhower.sph.harvard.edu/tools/) tools.
    - `run_KneadData.py`: Running [KneadData](https://huttenhower.sph.harvard.edu/kneaddata/) for read QC and host contaminants removal.
    - `run_MetaPhlAn4.py`: Running [MetaPhlAn4](https://huttenhower.sph.harvard.edu/metaphlan) for marker-gene based taxonomy profiling.
    - `run_HUMAnN3.py`: Running [HUMAnN3](https://huttenhower.sph.harvard.edu/humann) for functional profiling based on Uniref90 DB.
  - CAZyme: Script for carbohydrate-active enzyme (CAZyme) analysis.
    - `bbduk_adapter_trim.py`,  `bbduk_artifact_filter.py`, and `bbduk_quality_processing.py`: Running [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) for read QC.
    - `run_bowtie2.py` and `run_SAMtools.py`: Running [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [SAMtools](https://www.htslib.org/) for host-derived reads removal and contig coverage calculation.
    - `read_assembly.py`: Running [SPAdes](https://github.com/ablab/spades) for contig assembly.
    - `binning.py`: Binning metagenome-assembled genomes (MAGs) using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/), [MaxBin2](https://sourceforge.net/projects/maxbin2/), and [CONCOCT](https://github.com/BinPro/CONCOCT).
    - `bin_refinement.py`: Refining MAGs using [DAS Tool](https://github.com/cmks/DAS_Tool)
    - `binQC.py`: Assessing quality of MAGs using [CheckM](https://github.com/Ecogenomics/CheckM).
    - `cazyme_prediction.py`: Running [run_dbCAN](https://github.com/linnabrown/run_dbcan) for MAG CAZyme prediction.
    - `taxonomy_classification.py`: Running [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) for MAG taxonomy classification.
    - `gene_coverage_calculation.py`: Running `featureCounts` from [Subread](https://subread.sourceforge.net/) for CAZyme abundance calculation. The abundances are further normalized by gene length.
