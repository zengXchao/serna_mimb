#########################################################################
# Customizable Parameters
#########################################################################

# CPU_NUM: Sets the number of threads available for parallel computation
CPU_NUM=16

# ADAPT: Specifies the adapter sequence that needs to be trimmed from the 3' end of the FASTQ data
ADAPT="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

# CONV_SAMPLES: Specifies the sample names for conventional RNA-seq, with each line corresponding to a replicate. The sample names must end with '_repN' (where N=1, 2, ...).
CONV_SAMPLES="
	HEK_Conv_rep1
	HEK_Conv_rep2
	"

# CONV_FASTQ_FILES: Specifies the download URLs for the FASTQ files corresponding to the sample names listed in the CONV_SAMPLES parameter.
# If the data is paired-end, each line should specify the URLs for both Read1 and Read2 FASTQ files, separated by a comma. For example: R1_URL,R2_URL
CONV_FASTQ_FILES="
	https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA009/DRA009793/DRX205290/DRR214955.fastq.bz2
	https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA012/DRA012807/DRX310377/DRR321030.fastq.bz2
	"

# IMPR_SAMPLES: Specifies the sample names for improved (semi-extractable) RNA-seq, with each line corresponding to a replicate. The sample names must end with '_repN' (where N=1, 2, ...).
IMPR_SAMPLES="
	HEK_Impr_rep1 
	HEK_Impr_rep2
	"

# IMPR_FASTQ_FILES: Specifies the download URLs for the FASTQ files corresponding to the sample names listed in the IMPR_SAMPLES parameter.
# If the data is paired-end, each line should specify the URLs for both Read1 and Read2 FASTQ files, separated by a comma. For example: R1_URL,R2_URL
IMPR_FASTQ_FILES="
	https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA009/DRA009793/DRX205291/DRR214956.fastq.bz2
	https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA012/DRA012807/DRX310378/DRR321031.fastq.bz2
	"


#########################################################################
# Predefined Parameters
#########################################################################
DATA_DIR="01_data"
ALN_DIR="02_aln"
STRINGTIE_DIR="03_stringtie"
FEATURECOUNTS_DIR="04_featurecounts"
EDGER_DIR="05_edgeR"
CHMHMM_DIR="06_chromatinHMM"
APEX_DIR="06_apex"
RNA_DNA_DIR="08_RNA-DNA"
MFE_DIR="09_MFE"
LEN_DIR="10_len"
GC_DIR="11_GC"
TE_DIR="12_TE"
RNA_RNA_DIR="13_RNA-RNA"
MOTIF_DIR="14_motif"
SEEKR_DIR="15_seekr"
OUTPUT_DIR="output"

#########################################################################
# Data preparation
#########################################################################
mkdir -p 01_data/seq
CONV_SAMPLE_ARRAY=($CONV_SAMPLES)
CONV_FASTQ_ARRAY=($CONV_FASTQ_FILES)
for ((i=0; i<${#CONV_SAMPLE_ARRAY[@]}; i++))
do
	sample_name=${CONV_SAMPLE_ARRAY[$i]}
	sample_url=${CONV_FASTQ_ARRAY[$i]}
	if [[ $sample_url == *","* ]]; then
		IFS=',' read -r -a urls <<< "$sample_url"
		sample_url_r1=$(urls[0])
		sample_url_r2=$(urls[1])
		sample_file_r1=01_data/seq/${sample_name}_R1.fastq.bz2
		sample_file_r2=01_data/seq/${sample_name}_R2.fastq.bz2
		wget -O ${sample_file_r1} ${sample_url_r1}
		bzip2 -d ${sample_file_r1}
		wget -O ${sample_file_r2} ${sample_url_r2}
		bzip2 -d ${sample_file_r2}
		is_PE="TRUE"
	else
		sample_file=01_data/seq/${sample_name}.fastq.bz2
		wget -O ${sample_file} ${sample_url}
		bzip2 -d ${sample_file}
		is_PE="FALSE"
	fi	
done

IMPR_SAMPLE_ARRAY=($IMPR_SAMPLES)
IMPR_FASTQ_ARRAY=($IMPR_FASTQ_FILES)
for ((i=0; i<${#IMPR_SAMPLE_ARRAY[@]}; i++))
do
	sample_name=${IMPR_SAMPLE_ARRAY[$i]}
	sample_url=${IMPR_FASTQ_ARRAY[$i]}
	if [[ $sample_url == *","* ]]; then
		IFS=',' read -r -a urls <<< "$sample_url"
		sample_url_r1=$(urls[0])
		sample_url_r2=$(urls[1])
		sample_file_r1=01_data/seq/${sample_name}_R1.fastq.bz2
		sample_file_r2=01_data/seq/${sample_name}_R2.fastq.bz2
		wget -O ${sample_file_r1} ${sample_url_r1}
		bzip2 -d ${sample_file_r1}
		wget -O ${sample_file_r2} ${sample_url_r2}
		bzip2 -d ${sample_file_r2}
	else
		sample_file=01_data/seq/${sample_name}.fastq.bz2
		wget -O ${sample_file} ${sample_url}
		bzip2 -d ${sample_file}
	fi	
done


mkdir -p 01_data/ref
wget -O 01_data/ref/genome.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.p13.genome.fa.gz
gunzip 01_data/ref/genome.fa.gz
samtools faidx 01_data/ref/genome.fa
cut -f1,2 01_data/ref/genome.fa.fai > 01_data/ref/genome.seqsize
wget -O 01_data/ref/gencode.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip 01_data/ref/gencode.gtf.gz
wget --no-check-certificate -q -O 01_data/ref/rRNA.fa https://zengchao.me/files/human.rRNA.fa

mkdir -p 01_data/ref/starIndex
STAR --runMode genomeGenerate --genomeDir 01_data/ref/starIndex/gencode --genomeFastaFiles 01_data/ref/genome.fa --sjdbGTFfile 01_data/ref/gencode.gtf
STAR --runMode genomeGenerate --genomeDir 01_data/ref/starIndex/rRNA --genomeFastaFiles 01_data/ref/rRNA.fa
python ./bin/collapse_gencode_gtf.py --gtf 01_data/ref/gencode.gtf --out 01_data/ref/gencode.collapse.gtf
python ./bin/collapse_gencode_gtf.py --gtf 01_data/ref/gencode.gtf --includeIntron --out 01_data/ref/gencode.collapse.withintron.gtf
awk '{if($3!="transcript" && $3!="CDS"){printf "%s\tCOLLAPSE\ttranscript\t%s\t%s\t.\t%s\t.\tgene_id %s/%s transcript_id %s\n%s\tCOLLAPSE\texon\t%s\t%s\t.\t%s\t.\tgene_id %s/%s transcript_id %s\n",$1,$4,$5,$7,substr($10,1,length($10)-2),substr($12,2,length($12)),$14,$1,$4,$5,$7,substr($10,1,length($10)-2),substr($12,2,length($12)),$14}}' 01_data/ref/gencode.collapse.withintron.gtf > 01_data/ref/gencode.collapse.forstringtie.gtf


#########################################################################
# Detection of semi-extractable RNAs
#########################################################################

## Trimming adapter sequneces and mapping to rRNAs and the genome 
mkdir -p $ALN_DIR/01_trim $ALN_DIR/02_map2rrna $ALN_DIR/03_map2genome
for sample in $CONV_SAMPLES $IMPR_SAMPLES; do
	if [[ $is_PE == "FALSE" ]]; then
		cutadapt -j $CPU_NUM --overlap 5 --trim-n --max-n 1 --minimum-length 20 --quiet -a $ADAPT -o $ALN_DIR/01_trim/$sample.fastq $DATA_DIR/seq/$sample.fastq 
		STAR --genomeDir $DATA_DIR/ref/starIndex/rRNA --runThreadN $CPU_NUM \
			--outFilterMultimapNmax 1000 --outSAMmultNmax 1 --outReadsUnmapped Fastx \
			--readFilesIn $ALN_DIR/01_trim/$sample.fastq \
			--outFileNamePrefix $ALN_DIR/02_map2rrna/$sample.
		STAR --genomeDir $DATA_DIR/ref/starIndex/gencode --runThreadN $CPU_NUM --outFilterMultimapNmax 1 \
			--outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
			--readFilesIn $ALN_DIR/02_map2rrna/$sample.Unmapped.out.mate1 \
			--outFileNamePrefix $ALN_DIR/03_map2genome/$sample.
	else
		cutadapt -j $CPU_NUM --overlap 5 --trim-n --max-n 1 --minimum-length 50:50 --quiet -a $ADAPT -A $ADAPT \
		-o ${ALN_DIR}/01_trim/${sample}_R1.fastq -p ${ALN_DIR}/01_trim/${sample}_R2.fastq \
		${DATA_DIR}/seq/${sample}_R1.fastq ${DATA_DIR}/seq/${sample}_R2.fastq
		STAR --genomeDir $DATA_DIR/ref/starIndex/rRNA --runThreadN $CPU_NUM \
			--outFilterMultimapNmax 1000 --outSAMmultNmax 1 --outReadsUnmapped Fastx \
			--readFilesIn ${ALN_DIR}/01_trim/${sample}_R1.fastq ${ALN_DIR}/01_trim/${sample}_R2.fastq \
			--outFileNamePrefix $ALN_DIR/02_map2rrna/$sample.
		STAR --genomeDir $DATA_DIR/ref/starIndex/gencode --runThreadN $CPU_NUM --outFilterMultimapNmax 1 \
			--outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
			--readFilesIn $ALN_DIR/02_map2rrna/$sample.Unmapped.out.mate1 $ALN_DIR/02_map2rrna/$sample.Unmapped.out.mate2 \
			--outFileNamePrefix $ALN_DIR/03_map2genome/$sample.
	fi
done


## Transcriptome assembly
mkdir -p $STRINGTIE_DIR/tmp
for sample in $IMPR_SAMPLES; do
	stringtie -p $CPU_NUM --rf -o $STRINGTIE_DIR/tmp/$sample.gtf -l $sample $ALN_DIR/03_map2genome/$sample.Aligned.sortedByCoord.out.bam
	awk '{if($7=="+"){print $0}}' $STRINGTIE_DIR/tmp/$sample.gtf > $STRINGTIE_DIR/tmp/$sample.plus.gtf
	awk '{if($7=="-"){print $0}}' $STRINGTIE_DIR/tmp/$sample.gtf > $STRINGTIE_DIR/tmp/$sample.minus.gtf
	plus_gtfs="$plus_gtfs $STRINGTIE_DIR/tmp/$sample.plus.gtf"
	minus_gtfs="$minus_gtfs $STRINGTIE_DIR/tmp/$sample.minus.gtf"
	stringtie $ALN_DIR/03_map2genome/$sample.Aligned.sortedByCoord.out.bam -p $CPU_NUM -t -x --rf -e -G $DATA_DIR/ref/gencode.collapse.forstringtie.gtf -o $STRINGTIE_DIR/tmp/$sample.gtf
done
stringtie --merge -l STRG.plus -F 1 -o $STRINGTIE_DIR/tmp/plus.gtf $plus_gtfs
stringtie --merge -l STRG.minus -F 1 -o $STRINGTIE_DIR/tmp/minus.gtf $minus_gtfs
cat $STRINGTIE_DIR/tmp/plus.gtf $STRINGTIE_DIR/tmp/minus.gtf > $STRINGTIE_DIR/tmp/STRG.gtf
python bin/calculate_RS.py $STRINGTIE_DIR/tmp $IMPR_SAMPLES
python bin/merge_gtf.py $STRINGTIE_DIR/tmp/merged.withintron.gtf $DATA_DIR/ref/genome.seqsize $DATA_DIR/ref/gencode.collapse.gtf $DATA_DIR/ref/gencode.collapse.forstringtie.gtf $STRINGTIE_DIR/tmp/STRG.gtf $STRINGTIE_DIR/tmp/intron_rs.txt
awk '{if($2=="StringTie" && $3=="transcript"){printf "%s\t%d\t%s\t%s\t.\t%s\n",$1,$4-1,$5,substr($12,2,length($12)-3),$7}}' $STRINGTIE_DIR/tmp/merged.withintron.gtf| sort -k1,1 -k2,2n > $STRINGTIE_DIR/tmp/STRG.bed
bedtools cluster -s -d 10000 -i $STRINGTIE_DIR/tmp/STRG.bed > $STRINGTIE_DIR/tmp/STRG.cluster
awk '{if($2!="StringTie"){print}}' $STRINGTIE_DIR/tmp/merged.withintron.gtf > $STRINGTIE_DIR/tmp/merged.withoutSTR.gtf
python bin/cluster_STRG.py $STRINGTIE_DIR/tmp/STRG.cluster $STRINGTIE_DIR/tmp/merged.withoutSTR.gtf $STRINGTIE_DIR/merged.gtf
gffread -w $STRINGTIE_DIR/merged.fa -g $DATA_DIR/ref/genome.fa $STRINGTIE_DIR/merged.gtf
awk '$3=="exon" {print $1,$4,$5,$12":"$14,".",$7}' $STRINGTIE_DIR/merged.gtf | sed -e "s/\"//g" -e "s/;//g" -e "s/ /\t/g" | sort -k1,1 -k2,2n > $STRINGTIE_DIR/merged.bed
samtools faidx $STRINGTIE_DIR/merged.fa
rm -rf $STRINGTIE_DIR/tmp

## Quantification and differential expression analysis
mkdir -p $FEATURECOUNTS_DIR
for sample in $CONV_SAMPLES $IMPR_SAMPLES; do
	if [[ $is_PE == "FALSE" ]]; then
		featureCounts -T $CPU_NUM -g transcript_id -O -s 2 -a $STRINGTIE_DIR/merged.gtf -o $FEATURECOUNTS_DIR/$sample.txt $ALN_DIR/03_map2genome/$sample.Aligned.sortedByCoord.out.bam 2> /dev/null
	else
		featureCounts -T $CPU_NUM -g transcript_id -O -s 2 -p -a $STRINGTIE_DIR/merged.gtf -o $FEATURECOUNTS_DIR/$sample.txt $ALN_DIR/03_map2genome/$sample.Aligned.sortedByCoord.out.bam 2> /dev/null
	fi
done
mkdir -p $EDGER_DIR
CONV_LIST=$(echo "$CONV_SAMPLES" | tr -d '\n' | sed -e 's/[[:space:]]\+/ /g' | sed -e 's/ /,/g' -e 's/^,//' -e 's/,$//')
IMPR_LIST=$(echo "$IMPR_SAMPLES" | tr -d '\n' | sed -e 's/[[:space:]]\+/ /g' | sed -e 's/ /,/g' -e 's/^,//' -e 's/,$//')
python bin/runEdgeR.py -c $CONV_LIST -t $IMPR_LIST -f $FEATURECOUNTS_DIR -g $STRINGTIE_DIR/merged.gtf -o $EDGER_DIR

## Identification of semi-extractable RNAs and extractable RNAs
mkdir -p $OUTPUT_DIR
python bin/fig1.py $EDGER_DIR/diff_exp_tx.tsv 2 0.05 $STRINGTIE_DIR/merged.fa $STRINGTIE_DIR/merged.gtf $OUTPUT_DIR/fig1.png $OUTPUT_DIR/fig1_seRNAs.tsv $OUTPUT_DIR/fig1_exRNAs.tsv
awk '{print $1":"}' $OUTPUT_DIR/fig1_seRNAs.tsv | xargs -I x grep x $STRINGTIE_DIR/merged.bed | sort -k1,1 -k2,2n > $STRINGTIE_DIR/seRNA.bed
awk '{print $1":"}' $OUTPUT_DIR/fig1_exRNAs.tsv | xargs -I x grep x $STRINGTIE_DIR/merged.bed | sort -k1,1 -k2,2n > $STRINGTIE_DIR/exRNA.bed


#########################################################################
# Characterization of semi-extractable RNAs
#########################################################################

## Chromatin states
mkdir -p $CHMHMM_DIR
wget -O $CHMHMM_DIR/wgEncodeAwgSegmentationChromhmmHelas3.bed.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmHelas3.bed.gz
gunzip $CHMHMM_DIR/wgEncodeAwgSegmentationChromhmmHelas3.bed.gz
wget -O $CHMHMM_DIR/hg19ToHg38.over.chain.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
gunzip $CHMHMM_DIR/hg19ToHg38.over.chain.gz
liftOver $CHMHMM_DIR/wgEncodeAwgSegmentationChromhmmHelas3.bed $CHMHMM_DIR/hg19ToHg38.over.chain $CHMHMM_DIR/chromHMM.bed.tmp $CHMHMM_DIR/chromHMM.bed.unmapped
sort -k1,1 -k2,2n $CHMHMM_DIR/chromHMM.bed.tmp | cut -f1,2,3,4> $CHMHMM_DIR/chromHMM.bed
rm -rf $CHMHMM_DIR/wgEncodeAwgSegmentationChromhmmHelas3.bed $CHMHMM_DIR/hg19ToHg38.over.chain $CHMHMM_DIR/chromHMM.bed.tmp $CHMHMM_DIR/chromHMM.bed.unmapped
bedtools intersect -nonamecheck -wo -a $CHMHMM_DIR/chromHMM.bed -b $STRINGTIE_DIR/merged.bed > $CHMHMM_DIR/chromHMM_x_BG
bedtools intersect -nonamecheck -wo -a $CHMHMM_DIR/chromHMM.bed -b $STRINGTIE_DIR/seRNA.bed > $CHMHMM_DIR/chromHMM_x_SE
bedtools intersect -nonamecheck -wo -a $CHMHMM_DIR/chromHMM.bed -b $STRINGTIE_DIR/exRNA.bed > $CHMHMM_DIR/chromHMM_x_EX
python bin/table1.py $STRINGTIE_DIR/seRNA.bed $STRINGTIE_DIR/exRNA.bed $STRINGTIE_DIR/merged.bed $CHMHMM_DIR/chromHMM_x_SE $CHMHMM_DIR/chromHMM_x_EX $CHMHMM_DIR/chromHMM_x_BG $CHMHMM_DIR/summary.txt
cp $CHMHMM_DIR/summary.txt $OUTPUT_DIR/table1.tsv

## RNA subcellular location (APEX)
mkdir -p $APEX_DIR/data $APEX_DIR/aln $APEX_DIR/fpkm
fasterq-dump SRR7367268 -S -o $APEX_DIR/data/NU_target_rep1 
fasterq-dump SRR7367269 -S -o $APEX_DIR/data/NU_target_rep2 
fasterq-dump SRR7367270 -S -o $APEX_DIR/data/NU_control_rep1 
fasterq-dump SRR7367271 -S -o $APEX_DIR/data/NU_control_rep2 
fasterq-dump SRR7367272 -S -o $APEX_DIR/data/NC_target_rep1 
fasterq-dump SRR7367273 -S -o $APEX_DIR/data/NC_target_rep2 
fasterq-dump SRR7367274 -S -o $APEX_DIR/data/NC_control_rep1 
fasterq-dump SRR7367275 -S -o $APEX_DIR/data/NC_control_rep2 
fasterq-dump SRR7367276 -S -o $APEX_DIR/data/NL_target_rep1 
fasterq-dump SRR7367277 -S -o $APEX_DIR/data/NL_target_rep2 
fasterq-dump SRR7367278 -S -o $APEX_DIR/data/NL_control_rep1 
fasterq-dump SRR8646402 -S -o $APEX_DIR/data/NL_control_rep2 
fasterq-dump SRR7367280 -S -o $APEX_DIR/data/NP_target_rep1 
fasterq-dump SRR7367281 -S -o $APEX_DIR/data/NP_target_rep2 
fasterq-dump SRR7367282 -S -o $APEX_DIR/data/NP_control_rep1 
fasterq-dump SRR7367283 -S -o $APEX_DIR/data/NP_control_rep2 
fasterq-dump SRR7367284 -S -o $APEX_DIR/data/CY_target_rep1 
fasterq-dump SRR7367285 -S -o $APEX_DIR/data/CY_target_rep2 
fasterq-dump SRR7367286 -S -o $APEX_DIR/data/CY_control_rep1 
fasterq-dump SRR7367287 -S -o $APEX_DIR/data/CY_control_rep2 
fasterq-dump SRR7367288 -S -o $APEX_DIR/data/EM_target_rep1 
fasterq-dump SRR7367289 -S -o $APEX_DIR/data/EM_target_rep2 
fasterq-dump SRR7367290 -S -o $APEX_DIR/data/EM_control_rep1 
fasterq-dump SRR7367291 -S -o $APEX_DIR/data/EM_control_rep2 
fasterq-dump SRR7367296 -S -o $APEX_DIR/data/OM_target_rep1 
fasterq-dump SRR7367297 -S -o $APEX_DIR/data/OM_target_rep2 
fasterq-dump SRR7367298 -S -o $APEX_DIR/data/OM_control_rep1 
fasterq-dump SRR7367299 -S -o $APEX_DIR/data/OM_control_rep2 
fasterq-dump SRR7367300 -S -o $APEX_DIR/data/MM_target_rep1 
fasterq-dump SRR7367301 -S -o $APEX_DIR/data/MM_target_rep2 
fasterq-dump SRR7367303 -S -o $APEX_DIR/data/MM_control_rep1 
fasterq-dump SRR7367304 -S -o $APEX_DIR/data/MM_control_rep2 
fasterq-dump SRR7367292 -S -o $APEX_DIR/data/EL_target_rep1 
fasterq-dump SRR7367293 -S -o $APEX_DIR/data/EL_target_rep2 
fasterq-dump SRR7367294 -S -o $APEX_DIR/data/EL_control_rep1 
fasterq-dump SRR7367295 -S -o $APEX_DIR/data/EL_control_rep2 
for f in $APEX_DIR/data/*_1.fastq; do
	sample=`basename $f _1.fastq`
	STAR --genomeDir $DATA_DIR/ref/starIndex/gencode --runThreadN $CPU_NUM --outFilterMultimapNmax 1 \
		--outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
		--readFilesIn $APEX_DIR/data/${sample}_1.fastq $APEX_DIR/data/${sample}_2.fastq \
		--outFileNamePrefix $APEX_DIR/aln/$sample.
	stringtie -p $CPU_NUM --rf -e -G $STRINGTIE_DIR/merged.gtf -o $APEX_DIR/fpkm/$sample.gtf $APEX_DIR/aln/$sample.Aligned.sortedByCoord.out.bam
done
python bin/fig2.py $MFE_DIR/tmp/se_seq.tsv $MFE_DIR/tmp/ex_seq.tsv $APEX_DIR/fpkm $APEX_DIR/localization.png $APEX_DIR/localization_pval.tsv
cp $APEX_DIR/localization.png $OUTPUT_DIR/fig2.png
cp $APEX_DIR/localization_pval.tsv $OUTPUT_DIR/fig2_pval.tsv


## RNA-DNA interation (iMARGI)
mkdir -p $RNA_DNA_DIR
wget -O $RNA_DNA_DIR/iMARGI.bedpe.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3478nnn/GSM3478205/suppl/GSM3478205%5FHEK%5FiMARGI.bedpe.gz
gunzip $RNA_DNA_DIR/iMARGI.bedpe.gz

awk '$2>0&&$3>0&&$5>0&&$6>0 {print $1"\t"$2"\t"$3"\t"$7"_1\n"$4"\t"$5"\t"$6"\t"$7"_2"}' $RNA_DNA_DIR/iMARGI.bedpe | sort -k1,1 -k2,2n > $RNA_DNA_DIR/iMARGI.bed
bedtools intersect -nonamecheck -wo -b $RNA_DNA_DIR/iMARGI.bed -a $STRINGTIE_DIR/seRNA.bed > $RNA_DNA_DIR/SE_x_iMARGI
bedtools intersect -nonamecheck -wo -b $RNA_DNA_DIR/iMARGI.bed -a $STRINGTIE_DIR/exRNA.bed > $RNA_DNA_DIR/EX_x_iMARGI
bedtools intersect -nonamecheck -wo -b $RNA_DNA_DIR/iMARGI.bed -a $STRINGTIE_DIR/merged.bed > $RNA_DNA_DIR/BG_x_iMARGI
python bin/fig3.py $STRINGTIE_DIR/merged.fa.fai $STRINGTIE_DIR/seRNA.bed $STRINGTIE_DIR/exRNA.bed $STRINGTIE_DIR/merged.bed $RNA_DNA_DIR/SE_x_iMARGI $RNA_DNA_DIR/EX_x_iMARGI $RNA_DNA_DIR/BG_x_iMARGI $RNA_DNA_DIR/fig3.png $RNA_DNA_DIR/fig3.pval.txt 
cp $RNA_DNA_DIR/fig3.png $RNA_DNA_DIR/fig3.pval.txt $OUTPUT_DIR/

## MFE
mkdir -p $MFE_DIR/tmp
python bin/fig4.py $STRINGTIE_DIR/merged.fa $OUTPUT_DIR/fig1_seRNAs.tsv $OUTPUT_DIR/fig1_exRNAs.tsv $MFE_DIR/tmp $MFE_DIR/fig4.png $MFE_DIR/fig4_pval.txt
cp $MFE_DIR/fig4.png $MFE_DIR/fig4_pval.txt $OUTPUT_DIR/

## Length
mkdir -p $LEN_DIR
python bin/fig5.py $STRINGTIE_DIR/merged.fa.fai $OUTPUT_DIR/fig1_seRNAs.tsv $OUTPUT_DIR/fig1_exRNAs.tsv $LEN_DIR/fig5.png $LEN_DIR/fig5_pval.txt
cp $LEN_DIR/fig5.png $LEN_DIR/fig5_pval.txt $OUTPUT_DIR/

## GC content
mkdir -p $GC_DIR
python bin/fig6.py $MFE_DIR/tmp/all_seq.tsv $MFE_DIR/tmp/se_seq.tsv $MFE_DIR/tmp/ex_seq.tsv $GC_DIR/fig6.png $GC_DIR/fig6_pval.txt
cp $GC_DIR/fig6.png $GC_DIR/fig6_pval.txt $OUTPUT_DIR/

## Repeat
mkdir -p $TE_DIR
wget -O $TE_DIR/repeatmasker.hg38.lib20240131.out.gz https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz
gunzip $TE_DIR/repeatmasker.hg38.lib20240131.out.gz
awk 'NR>3 {print $5"\t"$6"\t"$7"\t"$11}' $TE_DIR/repeatmasker.hg38.lib20240131.out | sort -k1,1 -k2,2n > $TE_DIR/repeatmasker.bed
bedtools intersect -nonamecheck -wo -b $TE_DIR/repeatmasker.bed -a $STRINGTIE_DIR/seRNA.bed > $TE_DIR/SE_x_TE
bedtools intersect -nonamecheck -wo -b $TE_DIR/repeatmasker.bed -a $STRINGTIE_DIR/exRNA.bed > $TE_DIR/EX_x_TE
bedtools intersect -nonamecheck -wo -b $TE_DIR/repeatmasker.bed -a $STRINGTIE_DIR/merged.bed > $TE_DIR/BG_x_TE
python bin/fig7.py $STRINGTIE_DIR/merged.fa.fai $STRINGTIE_DIR/seRNA.bed $STRINGTIE_DIR/exRNA.bed $STRINGTIE_DIR/merged.bed $TE_DIR/SE_x_TE $TE_DIR/EX_x_TE $TE_DIR/BG_x_TE $TE_DIR/fig7.png $TE_DIR/fig7_pval.txt
cp $TE_DIR/fig7.png $TE_DIR/fig7_pval.txt $OUTPUT_DIR/

# RNA-RNA interation (RIC-seq and PARIS)
mkdir -p $RNA_RNA_DIR/tmp
wget -O $RNA_RNA_DIR/tmp/RICseq_hubRNA.xlsx https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2249-1/MediaObjects/41586_2020_2249_MOESM10_ESM.xlsx
in2csv $RNA_RNA_DIR/tmp/RICseq_hubRNA.xlsx --sheet "Supplementary Table 7" > $RNA_RNA_DIR/tmp/RICseq_hubRNA.csv
awk -F',' 'NR>3 {print $5}' $RNA_RNA_DIR/tmp/RICseq_hubRNA.csv | sort > $RNA_RNA_DIR/tmp/RICseq_hubRNA
grep "^>" $STRINGTIE_DIR/merged.fa | sed -e 's/>//' -e 's/_RI//' | sort | uniq > $RNA_RNA_DIR/tmp/all_tx
comm -12 $RNA_RNA_DIR/tmp/RICseq_hubRNA $RNA_RNA_DIR/tmp/all_tx > $RNA_RNA_DIR/RICseq_hubRNA.txt
wget -O $RNA_RNA_DIR/tmp/PARIS.txt.gz http://rise.life.tsinghua.edu.cn/download-data/rise_human_transcriptome.txt.gz
gunzip $RNA_RNA_DIR/tmp/PARIS.txt.gz
awk '$21=="HEK293T" {print $12"\n"$14}' $RNA_RNA_DIR/tmp/PARIS.txt | sort | uniq -c | awk '$1 > 20 {print $2}' > $RNA_RNA_DIR/tmp/PARIS_hubRNA
comm -12 $RNA_RNA_DIR/tmp/PARIS_hubRNA $RNA_RNA_DIR/tmp/all_tx > $RNA_RNA_DIR/PARIS_hubRNA.txt
cut -f1 $OUTPUT_DIR/fig1_seRNAs.tsv | sed 's/_RI//g' | sort > $RNA_RNA_DIR/seRNA.txt
cut -f1 $OUTPUT_DIR/fig1_exRNAs.tsv | sed 's/_RI//g' | sort > $RNA_RNA_DIR/exRNA.txt
python bin/fig8.py $RNA_RNA_DIR/seRNA.txt $RNA_RNA_DIR/RICseq_hubRNA.txt SE "hub RNAs (RIC-seq)" $MFE_DIR/tmp/all_seq.tsv $RNA_RNA_DIR/fig8A.png $RNA_RNA_DIR/fig8A_pval.txt
python bin/fig8.py $RNA_RNA_DIR/exRNA.txt $RNA_RNA_DIR/RICseq_hubRNA.txt EX "hub RNAs (RIC-seq)" $MFE_DIR/tmp/all_seq.tsv $RNA_RNA_DIR/fig8C.png $RNA_RNA_DIR/fig8C_pval.txt
python bin/fig8.py $RNA_RNA_DIR/seRNA.txt $RNA_RNA_DIR/PARIS_hubRNA.txt SE "hub RNAs (PARIS)" $MFE_DIR/tmp/all_seq.tsv $RNA_RNA_DIR/fig8B.png $RNA_RNA_DIR/fig8B_pval.txt
python bin/fig8.py $RNA_RNA_DIR/exRNA.txt $RNA_RNA_DIR/PARIS_hubRNA.txt EX "hub RNAs (PARIS)" $MFE_DIR/tmp/all_seq.tsv $RNA_RNA_DIR/fig8D.png $RNA_RNA_DIR/fig8D_pval.txt
cp $RNA_RNA_DIR/fig8?.png $RNA_RNA_DIR/fig8?_pval.txt  $OUTPUT_DIR/


# RBP-binding motif
mkdir -p $MOTIF_DIR/tmp
wget -O $MOTIF_DIR/tmp/motif_db.tgz https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.24.tgz
tar -xzf $MOTIF_DIR/tmp/motif_db.tgz -C $MOTIF_DIR/tmp
mv $MOTIF_DIR/tmp/motif_databases/CISBP-RNA/Homo_sapiens.meme $MOTIF_DIR/tmp
python bin/table2.py $MFE_DIR/tmp/se_seq.tsv $MFE_DIR/tmp/ex_seq.tsv $MOTIF_DIR/tmp/Homo_sapiens.meme $MOTIF_DIR/tmp
echo -e "RBP\tbin1\tbin2\tbin3\tbin4\tbin5\taverage" > $MOTIF_DIR/RBP_preference_full.tsv
awk '{\
    if (match($1, /\([^)]*\)/)) {\
        sub(/\(.*\)/, substr($1, RSTART+1, RLENGTH-2), $1)
    } \
    print\
}' $MOTIF_DIR/tmp/RBP_preference.txt | sort -k7,7nr >> $MOTIF_DIR/RBP_preference_full.tsv
head -n 11 $MOTIF_DIR/RBP_preference_full.tsv > $MOTIF_DIR/RBP_preference_top10_bottom10.tsv
tail -n 10 $MOTIF_DIR/RBP_preference_full.tsv >> $MOTIF_DIR/RBP_preference_top10_bottom10.tsv
cp $MOTIF_DIR/RBP_preference_top10_bottom10.tsv $OUTPUT_DIR/table2.tsv
cp $MOTIF_DIR/RBP_preference_full.tsv $OUTPUT_DIR/table2_full.tsv


# k-mer clustering
mkdir -p $SEEKR_DIR/tmp
cp $MOTIF_DIR/tmp/fimo/SE.fa $SEEKR_DIR/tmp/SE.fa
seekr_kmer_counts $SEEKR_DIR/tmp/SE.fa -o $SEEKR_DIR/tmp/SE.freq -k 6
seekr_pearson $SEEKR_DIR/tmp/SE.freq $SEEKR_DIR/tmp/SE.freq -o $SEEKR_DIR/tmp/SE.corr
python bin/table3.py $SEEKR_DIR/tmp/SE.corr $SEEKR_DIR/tmp/SE.cluster 0.5 $MFE_DIR/tmp/se_seq.tsv $MFE_DIR/tmp/ex_seq.tsv $MOTIF_DIR/tmp $SEEKR_DIR/tmp/SE_RBP_preference.tsv
sed "s/>//" $SEEKR_DIR/tmp/SE.cluster > $SEEKR_DIR/cluster.tsv
head -n 1 $SEEKR_DIR/tmp/SE_RBP_preference.tsv > $SEEKR_DIR/cluster_RBP_preference.tsv
awk 'NR>1 {\
    if (match($1, /\([^)]*\)/)) {\
        sub(/\(.*\)/, substr($1, RSTART+1, RLENGTH-2), $1)
    } \
    print\
}' $SEEKR_DIR/tmp/SE_RBP_preference.tsv | sort -k3,3nr >> $SEEKR_DIR/cluster_RBP_preference.tsv
























