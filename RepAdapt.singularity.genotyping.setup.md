This repo contains a code walkthrough for preparing containers and variables for genotyping using Singularity images for RepAdapt tools in the command line.
This code can be incorporated into bash scripts (or Nextflow or snakemake), or run directly, and operates under the assumption that you have:
1. Apptainer/Singularity installed
2. GNU parallel installed
3. Shell parameters for scripts specified either directly or in a parameters file. These  are:

```
projdir = project directory
projname = project name
genome = genome for species
genomedir =  location of genome assembly
species = species name
repadaptimages = location of Singularity images
bamfile = list of realigned.bam files to genotype
runname = genotyping run
```

The Singularity images required can be built with this command:

```
singularity pull fastp_0.20.1.sif https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0 

singularity pull bwa_0.7.17.sif https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8 

singularity pull samtools_1.16.1.sif https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0 

singularity pull gatk_4.4.1.sif https://depot.galaxyproject.org/singularity/gatk4:4.1.0.0--0

singularity pull picard_2.26.3.sif https://depot.galaxyproject.org/singularity/picard:2.26.3--hdfd78af_0 

singularity pull bedtools_2.27.1.sif https://depot.galaxyproject.org/singularity/bedtools:2.27.1--0 

singularity pull gatk_3.8.9.sif https://depot.galaxyproject.org/singularity/gatk:3.8--9 

singularity pull bcftools_1.16.sif https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1

singularity pull vcftools_1.16.sif https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321hdcf5f25_9  
```


