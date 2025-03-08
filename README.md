# ASFV-P72-Nanopore-Reference-Assembly

This repository includes the file P72_Nanopore_Script.py, a Python-based pipeline designed for the sequential concatenation and analysis of fastq.gz files produced by an Oxford Nanopore Technologies sequencer. The pipeline aligns the aggregated sequence files to an established ASFV P72 reference dataset. For the first 50 concatenated files, the pipeline determines the best reference (based on percent identity), sequence variation between the sample and the reference, coverage, and mapping quality metrics. This pipeline is a modification of the <a href = https://github.com/Global-ASFV-Research-Alliance/ASFV_Pipeline> ASFV de novo assembly pipeline </a> and accordingly, uses the same docker/singularity container.

## Dependencies

All software and library dependencies are contained within the docker/singularity image. See installation below. In addition to a folder containing your fastq.gz files, the working directory must contain the following files and folder.

### Files

1) MetaDataStorage.csv
2) MetadataNew.csv
3) P72_Nanopore_Script.py

### Folder

1) p72_References

    a) a curated set of ASFV P72 fasta files that have been indexed with bwa-mem2

## Installation

### Download the Docker Image

The Image can be downloaded and formatted for docker or singuilarity using the following commands:

#### Docker

First, open Docker Desktop. Then in Powershell enter the following command:

```powershell
docker pull garadock/asfv_denovo_assembly_pipeline:v06
```

#### Singularity

Open the terminal and enter the following command:

```terminal
singularity build --fakeroot ASFV_denovo_assembly_pipeline_v06.sif docker://garadock/asfv_denovo_assembly_pipeline:v06
```

## Instructions for use

### Set up MetadataNew.csv

MetadataNew.csv can be set up to run a single run or multiple consecutive runs. Please see the included MetadataNew.csv as an example.

### Running the Script (Docker)

Using the terminal, enter the directory that contains the files and folders detailed in dependency section.

#### Windows

1) Open docker desktop

2) Open powershell and paste the following commands:

    ```powershell
    $LinHomeMount = $HOME.substring(2).replace('\','/')
    $LinPWD = $pwd.path.substring(2).replace('\','/')
    docker run -it -v ${HOME}:$LinHomeMount -w $LinPWD garadock/asfv_denovo_assembly_pipeline:v06 /bin/bash
    ```

    ```powershell
    python3 P72_Nanopore_Script.py --metadata MetadataNew.csv
    ```

### Running the Script (Singularity)

example.sh script:

```bash
#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=string
#SBATCH --cpus-per-task=10
#SBATCH --job-name=CustomJobName
#SBATCH --time=UNLIMITED
#SBATCH --output=CustomLogName.log


module load singularity

File_DIR=/your/directory/containing/files/

export SINGULARITY_BIND=$File_DIR

singularity exec "/path/to/the/container/ASFV_denovo_assembly_pipeline_v06.sif" \
python3 P72_Nanopore_Script.py --metadata MetadataNew.csv
```

```shell
sbatch example.sh
```

## References

 1 ) <a href='https://www.mdpi.com/1999-4915/16/10/1522'> O'Donnell V, Spinard E, Xu L, Berninger A, Barrette RW, Gladue DP, Faburay B. Full-Length ASFV B646L Gene Sequencing by Nanopore Offers a Simple and Rapid Approach for Identifying ASFV Genotypes. Viruses. 2024 Sep 26;16(10):1522. doi: 10.3390/v16101522. PMID: 39459857; PMCID: PMC11512349. </a>

 2 ) <a href = https://pmc.ncbi.nlm.nih.gov/articles/PMC11359534/>Spinard E, Dinhobl M, Erdelyan CNG, O'Dwyer J, Fenster J, Birtley H, Tesler N, Calvelage S, Leijon M, Steinaa L, O'Donnell V, Blome S, Bastos A, Ramirez-Medina E, Lacasta A, Ståhl K, Qiu H, Nilubol D, Tennakoon C, Maesembe C, Faburay B, Ambagala A, Williams D, Ribeca P, Borca MV, Gladue DP. A Standardized Pipeline for Assembly and Annotation of African Swine Fever Virus Genome. Viruses. 2024 Aug 13;16(8):1293. doi: 10.3390/v16081293. PMID: 39205267; PMCID: PMC11359534. </a>
