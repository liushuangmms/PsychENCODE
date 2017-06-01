DIR=$1
PID=$2
CHR=$3

#required by farnam cluster
module load GCCcore/5.4.0 

cd ${DIR}

#subset vcf using expression id
vcf-subset -c select_brain_genotypeid.txt ${PID}_merged_${CHR}.dose.vcf.gz> ${PID}.${CHR}.sub_gen.vcf
bgzip ${PID}.${CHR}.sub_gen.vcf
tabix -p vcf ${PID}.${CHR}.sub_gen.vcf.gz

#filter the genotype
zcat ${PID}.${CHR}.sub_gen.vcf.gz | java -jar ~/software/snpEff/SnpSift.jar filter "((MAF>=0.01)&&(R2>=0.3))">${PID}.${CHR}.dose.MAF01.R6.vcf
bgzip ${PID}.${CHR}.dose.MAF01.R6.vcf
tabix -p vcf ${PID}.${CHR}.dose.MAF01.R6.vcf.gz
vcftools --gzvcf ${PID}.${CHR}.dose.MAF01.R6.vcf.gz --hwe 0.000001 --recode -c >${PID}.${CHR}.dose.MAF01.R6.HWEe6.vcf


#add chr to genotype and record header
awk '{ if($0 !~ /^#/) print "chr"$0; else if(match($0,/(##contig=<ID=)(.*)/,m)) print m[1]"chr"m[2]; else print $0 }' ${PID}.${CHR}.dose.MAF01.R6.HWEe6.vcf> ${PID}.${CHR}.eqtl.vcf
bcftools view -h ${PID}.${CHR}.sub_gen.vcf.gz> ${PID}.${CHR}.hdr.txt

#delete replicated variants then sort
sort -u -t$'\t' -k3,3 ${PID}.${CHR}.eqtl.vcf>${PID}.${CHR}.eqtl.u3.vcf
cat ${PID}.${CHR}.eqtl.u3.vcf | vcf-sort > ${PID}.${CHR}.eqtl.u3s.vcf

#replace the header
bcftools reheader -h ${PID}.${CHR}.hdr.txt ${PID}.${CHR}.eqtl.u3s.vcf -o ${PID}.${CHR}.eqtl.u3s.nh.vcf
bgzip ${PID}.${CHR}.eqtl.u3s.nh.vcf
tabix -p vcf ${PID}.${CHR}.eqtl.u3s.nh.vcf.gz
