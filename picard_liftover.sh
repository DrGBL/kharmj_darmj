cd /scratch/richards/guillaume.butler-laporte/WGS/DataFreeze3/

gzip -d -f grch37.vcf.gz
cat header.txt grch37.vcf | bgzip -c > grch37.for.lift.vcf.gz
tabix -p vcf grch37.for.lift.vcf.gz

java -jar /project/richards/guillaume.butler-laporte/bin/picard.jar LiftoverVcf \
    I=grch37.for.lift.vcf.gz \
    O=grch37.lifted_over.vcf \
    CHAIN=hg19ToHg38.over.chain \
    REJECT=grch37.rejected_variants.vcf \
    R=/scratch/richards/guillaume.butler-laporte/WGS/1000G.PCA/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    WARN_ON_MISSING_CONTIG=true \
    RECOVER_SWAPPED_REF_ALT=true

#now upload these two files back to own desktop to datafreeze3 folder common.variants subfolder
awk '!/##/ { print $1, $2, $3, $4, $5 }' grch37.lifted_over.vcf > grch37to38.maf.1.exome.lifted.over.txt
awk '!/##/ { print $1, $2, $3, $4, $5, $7, $8 }' grch37.rejected_variants.vcf > grch37to38.maf.1.exome.failed.lift.over.txt
