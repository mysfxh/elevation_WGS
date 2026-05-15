
###1

bcftools query -l /home/GYLW_WGS_data/03.variants/All.bgzf.vcf.gz

bcftools query -l /home/GYLW_WGS_data/15657-BN250716NJ01S02N2-BN250716NJ01S02N2/86samples.recode.vcf.gz 

bcftools view 03.variants/All.vcf.gz -O z -o 03.variants/All.bgzf.vcf.gz

bcftools reheader -s ../new_sample.txt 86samples.recode.vcf.gz -o 86samples.recode_renamed.vcf.gz


(WGS) root@server:/home/GYLW_WGS_data/15657-BN250716NJ01S02N2-BN250716NJ01S02N2# bcftools merge   --threads 200   -m none   -Oz   -o AllGYLW_merged.vcf.gz \
   /home/GYLW_WGS_data/15657-BN250716NJ01S02N2-BN250716NJ01S02N2/86samples.recode_renamed.vcf.gz\
  /home/GYLW_WGS_data/03.variants/All.bgzf.vcf.gz




# 1. 提取所有SNP（双等位基因）
nohup bcftools view -v snps -m2 -M2 /home/GYLW_WGS_data/15657-BN250716NJ01S02N2-BN250716NJ01S02N2/AllGYLW_merged.vcf.gz -Oz -o /home/xiongh/2026/Rana/Filter_data/AllGYLW.snp1.sitefiltered.vcf.gz &

nohup bash -c "vcftools \
--gzvcf /home/xiongh/2026/Rana/Filter_data/AllGYLW.snp1.sitefiltered.vcf.gz  \
--minDP 10 \
--max-missing 0.75 \
--maf 0.05 \
--min-alleles 2 \
--max-alleles 2 \
--recode \
--recode-INFO-all \
--stdout | bgzip -c > /home/xiongh/2026/Rana/Filter_data/AllGYLW_0.75.10filtered.vcf.gz" > /home/xiongh/2026/Rana/Filter_data/nvcftools_0.75.log 2>&1 & 


plink --vcf /home/xiongh/2026/Rana/Filter_data/AllGYLW_0.75.10filtered.vcf.gz \
  --allow-extra-chr \
  --double-id \
 --make-bed \
  --out /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10filtered


#written.
#20073669 variants loaded from .bim file.
#201 people (0 males, 0 females, 201 ambiguous) loaded from .fam.
#Ambiguous sex IDs written to
#/home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10filtered.nosex .
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 201 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.831572.
#20073669 variants and 201 people pass filters and QC.
#Note: No phenotypes present.
#--make-bed to




# 1. 重新设置 SNP ID，避免重复 ID
plink2 \
  --bfile /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10filtered \
  --set-all-var-ids @:#\$r,\$a \
  --new-id-max-allele-len 1000 truncate \
  --make-bed \
  --out /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10filtered.unique


# 2. 进行 LD pruning，生成 .prune.in 和 .prune.out
plink \
  --bfile /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10filtered.unique \
  --allow-extra-chr \
  --indep-pairwise 100 1 0.8 \
  --out /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10filtered.unique.100_1_0.8


# 3. 提取 LD 过滤后保留下来的 SNP
plink \
  --bfile /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10filtered.unique \
  --allow-extra-chr \
  --extract /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10filtered.unique.100_1_0.8.prune.in \
  --make-bed \
  --out /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8


# 4. 把 LD 过滤后的数据导出为 VCF
plink \
  --bfile /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8 \
  --allow-extra-chr \
  --recode vcf-iid bgz \
  --out /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked


cd /home/xiongh/2026/Rana/Filter_data/plink/

############################################################
# 1. 只保留 2–12 号染色体，去掉 1 号和所有 contig/scaffold
############################################################

awk 'BEGIN{OFS="\t"} 
  ($1 ~ /^[0-9]+$/ && $1 >= 2 && $1 <= 12) {print}
' AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.bim \
> AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.chr2to12.admix.bim


############################################################
# 2. 提取对应 SNP 的 bed/fam/bim
#    注意：不能只改 bim 然后直接 cp bed/fam
#    因为 bim 行数变少了，bed 也必须同步过滤
############################################################

awk '{print $2}' AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.chr2to12.admix.bim \
> keep.chr2to12.snps

plink \
  --bfile AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8 \
  --allow-extra-chr \
  --extract keep.chr2to12.snps \
  --make-bed \
  --out AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.chr2to12.admix


############################################################
# 3. 跑 ADMIXTURE
############################################################

cd /home/xiongh/2026/Rana/Filter_data/plink/

# 新建结果文件夹
mkdir -p /home/xiongh/2026/Rana/Filter_data/plink/admixture_chr2to12

# 跑 ADMIXTURE，并把结果存到新文件夹
nohup bash -c '
for K in $(seq 1 10); do
  admixture -j200 --cv AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.chr2to12.admix.bed $K \
  | tee /home/xiongh/2026/Rana/Filter_data/plink/admixture_chr2to12/AllGYLW_0.75.10.chr2to12.admixture.K${K}.log

  mv AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.chr2to12.admix.${K}.Q \
     /home/xiongh/2026/Rana/Filter_data/plink/admixture_chr2to12/AllGYLW_0.75.10.chr2to12.K${K}.Q

  mv AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.chr2to12.admix.${K}.P \
     /home/xiongh/2026/Rana/Filter_data/plink/admixture_chr2to12/AllGYLW_0.75.10.chr2to12.K${K}.P
done
' > /home/xiongh/2026/Rana/Filter_data/plink/admixture_chr2to12/adAllGYLW_0.75.10.chr2to12.log 2>&1 &

#######PCA 分析除去了性染色体的plink

plink \
  --bfile AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.chr2to12.admix \
  --allow-extra-chr \
  --pca 20 \
  --out ./PCa/AllGYLW_0.75.10.chr2to12.PCA


# 提取数据进行nj 树的构建
cd /home/xiongh/2026/Rana/Filter_data/plink/ 
touch NJ_data
cd NJ_data
#进行数据的提取
python ../vcf2phylip-master/vcf2phylip.py --input /home/xiongh/2026/Rana/Filter_data/plink/AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked.vcf.gz -f


#####第二种方法除去 性染色体的并合成vcf
cd /home/xiongh/2026/Rana/Filter_data/plink/

############################################################
# 1. 查看 VCF 里面的染色体编号
############################################################

bcftools index AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked.vcf.gz

bcftools query -f '%CHROM\n' AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked.vcf.gz \
  | sort -n | uniq -c \
  > chromosome.count.vcf.txt

cat chromosome.count.vcf.txt


############################################################
# 2. 去掉 1 号染色体，只保留 2–12 号染色体
############################################################

bcftools view \
  -r 2,3,4,5,6,7,8,9,10,11,12 \
  -Oz \
  -o AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked.chr2to12.vcf.gz \
  AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked.vcf.gz


############################################################
# 3. 为新的 VCF 建索引
############################################################

bcftools index AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked.chr2to12.vcf.gz


############################################################
# 4. 检查是否只剩 2–12 号染色体
############################################################

bcftools query -f '%CHROM\n' AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked.chr2to12.vcf.gz \
  | sort -n | uniq -c \
  > chromosome.count.vcf.chr2to12.txt

cat chromosome.count.vcf.chr2to12.txt


############################################################
# 5. 再用 vcf2phylip 转换为 fasta
############################################################

python ../vcf2phylip-master/vcf2phylip.py \
  --input AllGYLW_0.75.10.snapp.filtered.unique.100_1_0.8.unlinked.chr2to12.vcf.gz \
  -f
