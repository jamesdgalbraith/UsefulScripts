# An assorted collection of useful scripts for bioinformatic purposes

## Coverters: for converting weird formats to normal ones
trfToGff.py converts the -ngs output of Tandem Repeat Finder to GFF format. For example:

python trfToGff.py -i </path/to/trf/ngs/output> -o </path/to/out.gff>

irfRoBed.py converts the -ngs output of Inverted Repeat Finder to BED format. For example:
python irfToBed.py <irf_output> > irf.bed

## Calculators: for calculating various statistics
### Heterozygosity
het_windower.py counts the number heterozygous and homozygous sites over windows from a VCF called with mpileup froma single BAM (only used for single individual sequencing)

pooled_het_windower.py calculates the pooled heterozygosity over windows from a VCF called with mpileup froma single BAM (only used for pooled sequencing)

Before either of these steps variants should be using mpileup. For example:
```
# Align PacBio to genome, filter out secondary mappings
minimap2 -t ${THREADS} -ax map-pb ${GENOME} ${READS} | \
    samtools view -@ ${THREADS} -F 260 -b | \
    samtools sort -@ ${THREADS} > ${BASE}.bam
samtools index ${BASE}.bam
# Call variants
bcftools mpileup --threads ${THREADS} -Ou -f ${GENOME} ${BASE}.bam --max-depth 1000 --min-MQ 60 | \
	bcftools call --threads ${THREADS} -mA -V indels -Ob -o ${BASE}.mpileup.bcf
bcftools view -O v -V indels ${base}.mpileup.bcf --threads ${THREADS} -o ${BASE}.mpileup.vcf
```
