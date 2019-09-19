[fuzl@192.168.10.2 /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI
cp    /data/data/dataparsing/R100400180139/v300013477/L0*/UID_result_*/*HD*gz ./

for i in /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/*_1.fq.gz 
do
echo " perl  /home/fuzl/pipeline/DNA/Fastqc_bwa_gatk4_Mutect2_HC_freebayes_annovar_v3.2.pl -fq1 $i -fq2 ${i/_1/_2} \
-outdir ${i%%_1.fq.gz} -bed /home/fuzl/bed/162/162.slop150.sort.merge.bed  -gatktype 4  "
done >run.sh

nohup sh /home/fuzl/script/shell.sh run.sh 5 &

debug
单4的 reads ID中有"_"，比对出错

for i in *single*_1.fq.gz
do
echo "gunzip $i && sed -i 's/1_/1 /' ${i%%.gz} && gzip ${i%%.gz} &"
done 

for i in *single*_2.fq.gz
do
echo "gunzip $i && sed -i 's/2_/2 /' ${i%%.gz} && gzip ${i%%.gz} &"
done

for i in /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/*single*_1.fq.gz 
do
echo " perl  /home/fuzl/pipeline/DNA/Fastqc_bwa_gatk4_Mutect2_HC_freebayes_annovar_v3.2.pl -fq1 $i -fq2 ${i/_1/_2} \
-outdir ${i%%_1.fq.gz} -bed /home/fuzl/bed/162/162.slop150.sort.merge.bed  -gatktype 4  "
done >debug_run2.sh
#####################
#########################
#IDES
#########################
mkdir iDES_analysis
[fuzl@192.168.10.2 /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis 16:53 #839]
$ ln -s /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/*/GATK/*[e3].bam ./

mkdir output 
rename bam sorted.bam *bam 
for i in *.bam 
do
echo "nohup perl /home/fuzl/soft/iDES/ides-bam2freq.pl -o output -t 12 ./  /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta /home/fuzl/bed/162/162.slop150.sort.merge.bed &"
done >ides-bam2freq.sh

sh ides-bam2freq.sh
perl  /home/fuzl/soft/iDES/ides-makedb.pl   -o output/ -a capp-seq-db output/ -n 7 
for i in output/*Q30.txt ; do perl /home/fuzl/soft/iDES/ides-polishbg.pl -o output/ $i output/capp-seq-db_ides-bgdb.txt & 



#####################
去背景前后比较
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/output
for i in *Q30.rmbg.txt
do
#awk 'NR==FNR{a[$1"\t"$2]=$0}NR>FNR{if ($1"\t"$2 in a){print a[$1"\t"$2]"\t"$0}}' v300013477_L01-1-3HD753_single.sorted.freq.paired.Q30.rmbg.txt v300013477_L01-1-3HD753_single.sorted.freq.paired.Q30.txt
awk 'NR==FNR{a[$1"\t"$2]=$0}NR>FNR{if ($1"\t"$2 in a){print a[$1"\t"$2]"\t"$0}}' $i ${i%%.Q30.rmbg.txt}.Q30.txt > paste/${i%%.*}_Q30rmbg_vs_Q30.txt & 
done
输出：/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/output/paste
##################
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/mutation
建标准，
grep SNV v300013477_L01-1-3HD753_single.txt|cut -f 1-15>snv.contr

标准参照8个snv中加入降噪结果
for i in *txt 
do
awk 'NR==FNR{a[$1"\t"$2]=$0}NR>FNR{if (FNR==1){print "Chr\ts\te\tref\talt\tgene\tp\t"a["CHR\tPOS"]};{if ($1"\t"$2 in a){print $0"\t"a[$1"\t"$2]}}}'  /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/output/paste/${i%%.txt}_Q30rmbg_vs_Q30.txt snv.contr >snv_8hot/$i.Q30 
done

加入突变频率
for i in *txt
do
awk 'NR==FNR{a[$6"\t"$7]=$1"\t"$2"\t"$3}NR>FNR{if (FNR==1){print "Freq\tDP\tAO\t"$0}else{if ($1"\t"$2 in a){print a[$1"\t"$2]"\t"$0}else{print "-\t-\t-\t"$0}}}' $i snv_8hot/$i.Q30   >snv_8hot/$i.Q30.freq 
done

cd snv_8hot
grep ^ *freq |sed 's/:/\t/'  >total.cnv

统计indel
grep -v SNV v300013477_L01-1-3HD753_single.txt|cut -f 6-15 >indel/contr
cat indel/contr

for i in *txt 
do
awk 'NR==FNR{a[$6"\t"$7]=$1"\t"$2"\t"$3}NR>FNR{if (FNR==1){print "Freq\tDP\tAO"};{if ($1"\t"$2 in a){print a[$1"\t"$2]"\t"$0}else{print "-\t-\t-\t"$0}}}' $i indel/contr   >indel/$i.indel.freq 
done


#########################
#baseline 构建
##########################
519:
fuzl@192.168.10.2 /home/fuzl/soft/iDES/519MGI_baseline 
ln -s  /home/fuzl/project/huada_MGI_shenchan_liuzq/519/analysis/*B/GATK_pretreatment/*sort.bam ./
rename sort sorted *bam 
mkdir output 
#cp /home/fuzl/project/boke519_test9sample/508p.sort.merge.bed  /home/fuzl/bed/hanym_508p.sort.merge.bed

echo "nohup perl /home/fuzl/soft/iDES/ides-bam2freq.pl -o output -t 12 ./  /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta /home/fuzl/bed/hanym_508p.sort.merge.bed &"
 >ides-bam2freq.sh
sh ides-bam2freq.sh
##等待分析完成
perl  /home/fuzl/soft/iDES/ides-makedb.pl   -o output/ -a 519MGI-db output/ -n 19

###################
162:
cd /home/fuzl/soft/iDES/162_baseline 
mkdir output
rename mem.rmdup sort.bam *bam 
nohup perl /home/fuzl/soft/iDES/ides-bam2freq.pl -o output -t 12 ./  /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta /home/fuzl/bed/162/162.slop150.sort.merge.bed&
#for i in *.bam 
#do
echo "nohup perl /home/fuzl/soft/iDES/ides-bam2freq.pl -o output -t 12 ./  /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta /home/fuzl/bed/162/162.slop150.sort.merge.bed &" >
>ides-bam2freq.sh
sh ides-bam2freq.sh
perl  /home/fuzl/soft/iDES/ides-makedb.pl   -o output/ -a 162-db output/ -n 15

############
#cappseq 143 
/home/fuzl/soft/iDES/143_baseline  #
echo "nohup perl /home/fuzl/soft/iDES/ides-bam2freq.pl -o output -t 12 ./  /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta /home/fuzl/bed/143gene.bed.bed &"
>ides-bam2freq.sh
perl  /home/fuzl/soft/iDES/ides-makedb.pl   -o output/ -a 143-db output/ -n 21

###################################################

162 MGI分析
mv output input
for i in /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/input/*Q30.txt 
do 
perl /home/fuzl/soft/iDES/ides-polishbg.pl -o /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/output_162baseline/ $i /home/fuzl/soft/iDES/162_baseline/output/162-db_ides-bgdb.txt  & 
done

去背景前后比较
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/output_162baseline/
mkdir paste
for i in *Q30.rmbg.txt
do
#awk 'NR==FNR{a[$1"\t"$2]=$0}NR>FNR{if ($1"\t"$2 in a){print a[$1"\t"$2]"\t"$0}}' v300013477_L01-1-3HD753_single.sorted.freq.paired.Q30.rmbg.txt v300013477_L01-1-3HD753_single.sorted.freq.paired.Q30.txt
awk 'NR==FNR{a[$1"\t"$2]=$0}NR>FNR{if ($1"\t"$2 in a){print a[$1"\t"$2]"\t"$0}}' $i ../input/${i%%.Q30.rmbg.txt}.Q30.txt > paste/${i%%.*}_Q30rmbg_vs_Q30.txt & 
done
输出：/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/output_162baseline/paste

标准参照8个snv中加入162基线降噪结果
cd /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/mutation
for i in *txt 
do
awk 'NR==FNR{a[$1"\t"$2]=$0}NR>FNR{if (FNR==1){print "Chr\ts\te\tref\talt\tgene\tp\t"a["CHR\tPOS"]};{if ($1"\t"$2 in a){print $0"\t"a[$1"\t"$2]}}}' \
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/output_162baseline/paste/${i%%.txt}_Q30rmbg_vs_Q30.txt \
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/mutation/snv.contr >/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/mutation/snv_8hot_162baseline/$i.Q30 
done

加入突变频率
for i in *txt
do
awk 'NR==FNR{a[$6"\t"$7]=$1"\t"$2"\t"$3}NR>FNR{if (FNR==1){print "Freq\tDP\tAO\t"$0}else{if ($1"\t"$2 in a){print a[$1"\t"$2]"\t"$0}else{print "-\t-\t-\t"$0}}}' \
$i snv_8hot_162baseline/$i.Q30  >snv_8hot_162baseline/$i.Q30.freq 
done




samtools view v300013477_L02-1-3HD753/GATK/v300013477_L02-1-3HD753.bam -q 30 -L mutation/snv.contr.bed |grep -v MD:Z:100  >mutation/snv_8hot_mismash.sam
UID 去重：#*_1.fq.gz

for i in *_1.fq.gz ## v300013477_L01-1-3HD753_single_1.fq.gz 
#for i in  v300013477_L02-1-3HD753_1.fq.gz 
do
awk 'NR==FNR{if (NR%4==1){a[$1]=$NF}}NR>FNR{if (/^@/){print $0}else{if ("@"$1"/1" in a){print $0"\tID:Z:"a["@"$1"/1"]}}}' <(zcat $i ) <(samtools view -h ${i%%_1.fq.gz}/GATK/${i%%_1.fq.gz}.bam) \
|awk  '{if($NF == CUST_ID ){$2+=1024;print $0 }else{if($NF != CUST_ID){CUST_ID=$NF;print $0}}}' OFS="\t" |samtools view -Sb >${i%%_1.fq.gz}/GATK/${i%%_1.fq.gz}.uid_make.bam &
done

awk 'NR==FNR{if (NR%4==1){a[$1]=$NF}}NR>FNR{if (/^@/){print $0}else{if ("@"$1"/1" in a){print $0"\tID:Z:"a["@"$1"/1"]}}}' <(zcat /data/fuzl/project/huada_MGI_shenchan_liuzq_UMI/v300013477_L01-1-3HD753/trim/v300013477_L01-1-3HD753_1_paired.fq.gz ) <(samtools view -h /data/fuzl/project/huada_MGI_shenchan_liuzq_UMI/v300013477_L01-1-3HD753/GATK/v300013477_L01-1-3HD753.bam)|awk  '{if($NF == CUST_ID ){$2+=1024;print $0 }else{if($NF != CUST_ID){CUST_ID=$NF;print $0}}}' OFS="\t" |samtools view -Sb >/data/fuzl/project/huada_MGI_shenchan_liuzq_UMI/v300013477_L01-1-3HD753/GATK/v300013477_L01-1-3HD753.uid_make.bam 
#攀的生产 python /home/lip/0626/reads_bases.py /home/lip/0626/upload_13.txt /data/data/task/ /home/lip/0626/vcf/ /home/lip/0626/result/ /data/data/backup/history/ /report/2019-8-22数据统计nextseq-190133/
==========

for i in /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/*_1.fq.gz 
do
echo " perl  /home/fuzl/pipeline/DNA/Fastqc_bwa_gatk4_Mutect2_HC_freebayes_annovar_v4.pl -fq1 $i -fq2 ${i/_1/_2} \
-outdir ${i%%_1.fq.gz} -bed /home/fuzl/bed/162/162.slop150.sort.merge.bed  -gatktype 4  -uid "
done >run_uid.sh


for i in mutation/snv_8hot_162baseline/v300013477_L02-1-3HD753.txt.Q30.freq # mutation/snv_8hot_162baseline/*txt.Q30.freq
do
a=${i##*/}
b=${a%%.txt*}
awk 'NR==FNR{a[$6"\t"$7"\t"$9"\t"$10]=$1"\t"$2"\t"$3}NR>FNR{if (FNR==1){print "uid_Freq\tuid_DP\tuid_AO\tFreq\tDP\tAO\t"$0}else{if ($4"\t"$5"\t"$7"\t"$8 in a){print a[$4"\t"$5"\t"$7"\t"$8]"\t"$0}else{print "-\t-\t-\t"$0}}}' \
$b/uid/vcf/$b.hg19_multianno.txt.freq $i >mutation/snv_8hot_162baseline/$b.uid_vs_picard.Q30.freq 
done
cd /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid 
mkdir output 
rename _make .sorted *bam 
echo "nohup perl /home/fuzl/soft/iDES/ides-bam2freq.pl -o output -t 12 ./  /home/fuzl/soft/GATK/resources/bundle/hg19/ucsc.hg19.fasta /home/fuzl/bed/162/162.slop150.sort.merge.bed &" >
>ides-bam2freq.sh

用519panel降噪

for i in /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output/*Q30.txt 
do 
perl /home/fuzl/soft/iDES/ides-polishbg.pl -o /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output_519baseline/ $i \
/home/fuzl/soft/iDES/519MGI_baseline/output/519MGI-db_ides-bgdb.txt  & 
done

计算大于cutoff的突变
for i in /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output/*Q30.txt /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output_519baseline/*Q30.rmbg.txt
do
perl /home/fuzl/soft/iDES/ides-freq2vcf.pl -input $i -output $i.cutoff0.01_freq -cutoff 0.01 &
done

[fuzl@192.168.10.2 /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid 18:59 #253]
cp  ../mutation/snv_8hot_162baseline/* 8snv_control/

/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/output_519baseline



#####################
#519baseline+UID 去背景前后比较
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/
mkdir befor_vs_after_519base_cutoff0.001
for i in output/*cutoff0.001_freq
do
a=${i##*/}
awk 'NR==FNR{a[$1"\t"$2"\t"$15]=$0}NR>FNR{if ($1"\t"$2"\t"$15 in a){print a[$1"\t"$2"\t"$15]"\t"$0}}' $i output_519baseline/${a%%.txt*}.rmbg.txt.cutoff0.001_freq\
> befor_vs_after_519base_cutoff0.001/${a%%.*}_Q30_vs_Q30rmbg.freq & 
done
#输出：/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/befor_vs_after_519base_cutoff0.001

#标准参照8个snv中加入519基线降噪结果
cd /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/befor_vs_after_519base_cutoff0.001
for i in *Q30_vs_Q30rmbg.freq 
do
awk 'NR==FNR{a[$1"\t"$2"\t"$15]=$0}NR>FNR{if (FNR==1){print "Chr\ts\te\tref\talt\tgene\tp\t"a["CHR\tPOS\tALT"]};{if ($1"\t"$2"\t"$5 in a){print $0"\t"a[$1"\t"$2"\t"$5]}}}' \
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/befor_vs_after_519base_cutoff0.001/$i \
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/mutation/snv.contr >/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/8snv_control_519baseline/$i.8snv_control
done

加入突变频率
for i in *txt
do
awk 'NR==FNR{a[$6"\t"$7]=$1"\t"$2"\t"$3}NR>FNR{if (FNR==1){print "Freq\tDP\tAO\t"$0}else{if ($1"\t"$2 in a){print a[$1"\t"$2]"\t"$0}else{print "-\t-\t-\t"$0}}}' \
$i snv_8hot_162baseline/$i.Q30  >snv_8hot_162baseline/$i.Q30.freq 
done







for i in /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/output_162baseline/*.txt /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/input/*Q30.txt
do
perl /home/fuzl/soft/iDES/ides-freq2vcf.pl -input $i -output $i.cutoff0.001_freq -cutoff 0.001 &
done
mkdir befor_vs_after_162base_cutoff0.001
for i in input/*cutoff0.001_freq
do
a=${i##*/}
awk 'NR==FNR{a[$1"\t"$2"\t"$15]=$0}NR>FNR{if ($1"\t"$2"\t"$15 in a){print a[$1"\t"$2"\t"$15]"\t"$0}}' $i output_162baseline/${a%%.txt*}.rmbg.txt.cutoff0.001_freq\
> befor_vs_after_162base_cutoff0.001/${a%%.*}_Q30_vs_Q30rmbg.freq & 
done
#输出：/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis_uid/befor_vs_after_162base_cutoff0.001

#标准参照8个snv中加入519基线降噪结果
cd /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/befor_vs_after_162base_cutoff0.001
for i in *Q30_vs_Q30rmbg.freq 
do
awk 'NR==FNR{a[$1"\t"$2"\t"$15]=$0}NR>FNR{if (FNR==1){print "Chr\ts\te\tref\talt\tgene\tp\t"a["CHR\tPOS\tALT"]};{if ($1"\t"$2"\t"$5 in a){print $0"\t"a[$1"\t"$2"\t"$5]}}}' \
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/befor_vs_after_162base_cutoff0.001/$i \
/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/mutation/snv.contr >/home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/iDES_analysis/8snv_control_162baseline/$i.8snv_control
done



将热点的vcf挑出来，找到该点的base和uid
perl /home/fuzl/soft/iDES/ides-location-cnv-uid.pl -input v300013477_L04-1-24HD753/GATK/v300013477_L04-1-24HD753.uid_make.bam -hotspot_vcf mutation/snv.contr -output test_v2

#awk '{print $3"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}'  test|sort -k 2  |less 加到/home/fuzl/soft/iDES/ides-location-cnv-uid.pl了
for i in  */GATK/*uid_make.bam
do
perl /home/fuzl/soft/iDES/ides-location-cnv-uid.pl -input $i -hotspot_vcf mutation/snv.contr -output uid_hotspot/${i%%/*} & 
done 




test 
cd /home/fuzl/project/huada_MGI_shenchan_liuzq_UMI/test_pipeline_DNA_UID_v4
perl /home/fuzl/pipeline/DNA/Fastqc_bwa_gatk4_Mutect2_HC_freebayes_annovar_v4.pl -fq1 L02-1-3HD753_1.fq.gz -fq2 L02-1-3HD753_2.fq.gz -outdir L02-1-3HD753  -gatktype 4 -bed /home/fuzl/bed/162/162.slop150.sort.merge.bed -fusion -uid
perl /home/fuzl/pipeline/DNA/Fastqc_bwa_gatk4_Mutect2_HC_freebayes_annovar_v4.pl -fq1 L01-1-3HD753_single_1.fq.gz -fq2 L01-1-3HD753_single_2.fq.gz -outdir L01-1-3HD753_single -gatktype 4 -bed /home/fuzl/bed/162/162.slop150.sort.merge.bed -fusion -uid -single &


L02   110M  663s
L01   2.5G  46121s