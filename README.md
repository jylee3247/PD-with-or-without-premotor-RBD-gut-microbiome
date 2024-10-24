# Introduction
This repository includes custom codes used for the analyses of this manuscripts.<br/>
The codes have been tested on Ubuntu20.04
- Title: Distinct gut microbiome characteristics in patients with Parkinson’s disease based on the presence of premotor rapid-eye movement sleep behavior disorders
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
  - Metagenome-assembled genome (MAG) reconstruction: Script for MAG analysis.
    - `bbduk_adapter_trim.py`,  `bbduk_artifact_filter.py`, and `bbduk_quality_processing.py`: Running [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) for read QC.
    - `run_bowtie2.py` and `run_SAMtools.py`: Running [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [SAMtools](https://www.htslib.org/) for host-derived reads removal and contig coverage calculation.
    - `read_assembly.py`: Running [SPAdes](https://github.com/ablab/spades) for contig assembly.
    - `binning.py`: Binning metagenome-assembled genomes (MAGs) using [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/), [MaxBin2](https://sourceforge.net/projects/maxbin2/), and [CONCOCT](https://github.com/BinPro/CONCOCT).
    - `bin_refinement.py`: Refining MAGs using [DAS Tool](https://github.com/cmks/DAS_Tool)
    - `binQC.py`: Assessing quality of MAGs using [CheckM](https://github.com/Ecogenomics/CheckM).
    - `taxonomy_classification.py`: Running [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) for MAG taxonomy classification.
    - `gene_coverage_calculation.py`: Running `featureCounts` from [Subread](https://subread.sourceforge.net/) for CAZyme abundance calculation. The abundances are further normalized by gene length.
  - CAZyme: Script for carbohydrate-active enzyme (CAZyme) analysis.
    - `cazyme_prediction.py`: Running [run_dbCAN](https://github.com/linnabrown/run_dbcan) for MAG CAZyme prediction.
  - `MAG_CsgA_with_reference_seqs.msa` : A multiple sequence alignment (MSA) file with fasta format. CsgA amino acid sequences were retrevied from the MAGs obatined from current study, and aligned by [Clastal Omega](http://www.clustal.org/omega/) with reference CsgA sequences (UniProt accession numbers `P28307` for *Escherichia coli* and `A0A9Q7ZLG0` for *Citrobacter youngae*).

* **PD-gut-microbiome-dynamics**: Analysis for exploring gut microbiome dynamics of patients with PD according to presence of the premotor RBD and disease progression.
  - `PD-gut-microbiome-dynamics.R`: R script used for the analysis.
  - `16S.RData`: 16S rRNA gene amplicon sequencing data for the analysis.
  - `WGS.RData`: whole shotgun metagenome sequencing data for the analysis.
  - `utilities.RData`: Utility codes for the analysis.

*All the raw data from the differential abundance analyses are available from the `Additional file 3`.*
