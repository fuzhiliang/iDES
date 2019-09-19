MGI_UID_split
/data/data/dataparsing/V300013665

#L03：
cd /data/data/dataparsing/R100400180139/v300013477/L03

#perl /data/data/dataparsing/V300013665/SplitDualBarcodes_fuzl.pl -r1 v300013477_L02_read_1.fq.gz -r2 v300013477_L02_read_2.fq.gz -e 0 -f 101 -b uid.txt

cat >L03_reverse_barcord.txt
GACGAACT        ACTCGATC
TACTGCTC        TGAGCTGT
TACTGCTC        TGAGCTGT
GACGAACT        ACTCGATC
lip@vm1:/data/data/dataparsing/R100400180139/v300013477/L03$ sh  UID_result_v3/*sh *



双8
perl /home/lip/fuzl/bin/MGI_SplitDualBarcodes.pl -r1 v300013477_L03_read_1.fq.gz  -r2 v300013477_L03_read_2.fq.gz  -e 0  -f 101 -b L03_reverse_barcord.txt  -o  UID_result_v3/ 

单4+4

perl /home/lip/fuzl/bin/MGI_SplitBarcode_uid.pl -r1 v300013477_L03_read_1.fq.gz -r2 v300013477_L03_read_2.fq.gz -f 109 -l 112 -b L03_4bp_barcord.txt -o UID_result_single/ &



L01
perl /home/lip/fuzl/bin/MGI_SplitBarcode_uid.pl -r1 v300013477_L01_read_1.fq.gz -r2 v300013477_L01_read_2.fq.gz -f 109 -l 112 -b L01_4bp_barcond.txt -o UID_result_single/ &

L02
perl /home/lip/fuzl/bin/MGI_SplitDualBarcodes.pl -r1 v300013477_L02_read_1.fq.gz  -r2 v300013477_L02_read_2.fq.gz   -f 101 -b L02_barcord.txt  -o  UID_result_v3/ 

L04
perl /home/lip/fuzl/bin/MGI_SplitDualBarcodes.pl -r1 v300013477_L04_read_1.fq.gz  -r2 v300013477_L04_read_2.fq.gz   -f 101 -b L04_barcord.txt  -o  UID_result_v3/ 



CAGTAACCCCAAAAACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCACTTTTAACAGTCACCCCCCAACTAACA




