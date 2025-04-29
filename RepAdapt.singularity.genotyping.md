This repo contains the code required for running genotyping steps with Singularity images. Setup instructions are [here](https://github.com/RepAdapt/singularity/blob/main/RepAdapt.singularity.genotyping.setup.md)

Trim with fastp:
```
#export variables to make parallel happy
export fastp
export projdir
export repadaptimages
export set

cat $projdir/sets/$set | \
  parallel -j 16 ' apptainer exec -B $projdir/reads \
  -B $projdir/trim \
  -B $projdir \
  $repadaptimages/fastp_0.20.1.sif  \
  bash -c " fastp -i $projdir/reads/{}_R1.fastq.gz \
  -I $projdir/reads/{}_R2.fastq.gz \
  -o $projdir/trim/{}_R1.trimmed.fastq.gz \
  -O $projdir/trim/{}_R2.trimmed.fastq.gz \
  --thread 4 "'
```

Indexing genome, samtools faidx, and sequence dictionary:

```
apptainer exec -B $genomedir \
$repadaptimages/bwa_0.7.17.sif bash -c "bwa index $genome"
apptainer exec -B $genomedir \
$repadaptimages/samtools_1.16.1.sif bash -c "samtools faidx $genome"
apptainer exec -B $genomedir \
$repadaptimages/gatk_4.4.1.sif bash -c "gatk CreateSequenceDictionary -R $genome"
```

