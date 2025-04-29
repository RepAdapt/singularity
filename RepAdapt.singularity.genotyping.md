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

align:

```

RGID=$(echo $ind |  sed 's/i5.*/i5/') ;
SMID=$(echo $ind | sed 's/NS.*i5.//') ;
LBID=$(echo $ind | sed 's/.UDP.*//');

apptainer exec \
  -B $projdir \
  -B $projdir/trim \
  -B $projdir/align \
  -B $genomedir \
  $repadaptimages/bwa_0.7.17.sif bash -c "bwa mem \
  -t 64 \
  -R \"@RG\tID:$RGID\tSM:$SMID\tLB:$LBID\" \
  $genome \
  ${ind}_R1.trimmed.fastq.gz ${ind}_R2.trimmed.fastq.gz" | \
   apptainer exec\
   -B $projdir \
   -B $projdir/align \
   -B $projdir/trim \
   $repadaptimages/samtools_1.16.1.sif bash -c " samtools sort \
   -o $projdir/align/${ind}.sorted.bam -T $ind \
   -@ 64 -m 2G"
```

Deduplicate aligned reads and index:

```

export projdir
export repadaptimages

#deduplicate in parallel

cat $projdir/sets/$set | \
  parallel --tmpdir $projdir/gatktemp \
  --jobs 14 \
  ' JAVA_OPTS="-Xmx5g" \
  apptainer exec -B $projdir/align \
   $repadaptimages/picard_2.26.3.sif   \
  bash -c "picard MarkDuplicates \
  INPUT=align/{}.sorted.bam \
  OUTPUT=align/{}.deDup.bam M=align/{}_deDupMetrics.txt \
  REMOVE_DUPLICATES=true \
  TMP_DIR=$projdir/gatktemp "'


cat $projdir/sets/$set | \
   parallel --jobs 25 \
   'apptainer exec\
   -B $projdir \
   -B $projdir/align \
   $repadaptimages/samtools_1.16.1.sif bash -c " samtools index align/{}.deDup.bam "'

```

Indel realignment

```
export repadaptimages
export projdir
export genome
export genomedir

cat $projdir/sets/$set | \
   parallel --tmpdir $projdir/gatktemp \
   --jobs 64 \
  ' JAVA_OPTS="-Xmx5g -Djava.io.tmpdir=/mnt/gatktemp" \
  apptainer exec \
   -B $projdir/gatktemp \
   -B $projdir/align \
   -B $genomedir \
   -B $projdir/gatktemp:/tmp \
   $repadaptimages/gatk_3.8.9.sif   \
  bash -c " gatk3  \
  -T RealignerTargetCreator \
  -R $genome \
  -I align/{}.deDup.bam \
  -o align/{}.intervals " '

cat $projdir/sets/$set | \
   parallel --tmpdir $projdir/gatktemp \
   --jobs 64 \
  ' JAVA_OPTS="-Xmx5g -Djava.io.tmpdir=/mnt/gatktemp" \
  apptainer exec \
   -B $projdir/gatktemp \
   -B $projdir/align \
   -B $genomedir \
   -B $projdir/gatktemp:/tmp \
   $repadaptimages/gatk_3.8.9.sif   \
  bash -c " gatk3  \
  -T IndelRealigner \
  -R $genome \
  -I align/{}.deDup.bam \
  -targetIntervals align/{}.intervals \
  -o align/{}.realigned.bam " '
```

Genotype by chromosome:

```
singularity exec \
  -B $projdir \
  -B $genomedir \
  $repadaptimages/bcftools_1.16.sif \
  bash -c " bcftools mpileup \
  -Ou \
  --threads 48  \
  -f $genome \
  --bam-list $bamfile -q 5 \
  -r $chrom \
  -a DP,AD \
 -I | \
  bcftools call \
  --threads 48 \
  -mv \
  -Ob \
  -f GQ,GP \
  -o $projdir/bcftools_out/$species.$projname.$runname.$chrom.bcf
"
```




