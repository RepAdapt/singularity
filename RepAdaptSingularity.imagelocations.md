# Singularity/Apptainer image locations and run commands
Activate singularity with:

```
mamba activate nf-core
```
through nf-core or :

```
mamba activate apptainer
```


## fastp (v.0.20.1)
```
singularity run https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0 
fastp
```
## bwa-mem (v.0.7.17-r1188)
```
singularity run https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8 
bwa index
bwa mem
```

## samtools (v.1.16.1) 
```
singularity run https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0 
samtools
```
## Picard Tools (v.2.26.3) 
```
singularity run https://depot.galaxyproject.org/singularity/picard:2.26.3--hdfd78af_0 
picard
```
## GATK v.3.8
```
singularity run https://depot.galaxyproject.org/singularity/gatk:3.8--9 
gatk3
```
## bcftools v. 1.16 
```
singularity run https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1
bcftools
```
## vcftools 0.1.16
```
singularity run https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321hdcf5f25_9  
vcftools
```
## BEDtools (v.2.27.1)
```
singularity run https://depot.galaxyproject.org/singularity/bedtools:2.27.1--0 
bedtools
```
## VarScan (v.2.4.2) 
```
singularity run https://depot.galaxyproject.org/singularity/varscan:2.4.2--0 
varscan
```
## GATK (v.4.1.0.0)
```
singularity run https://depot.galaxyproject.org/singularity/gatk4:4.1.0.0--0
gatk
```
