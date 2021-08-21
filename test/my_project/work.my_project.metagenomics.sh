ln -sb /USB/meta_genomics_pipeline/test/test_data/LV20R4_subsample.fastq.gz /USB/meta_genomics_pipeline/test/my_project/fastqc_result/sample01.raw.fastq.gz 
echo "finshed link fastq files" 
kneaddata --input /USB/meta_genomics_pipeline/test/my_project/fastqc_result/sample01.raw.fastq.gz --reference-db /USB/meta_genomics_pipeline/database/kneaddata_database --output /USB/meta_genomics_pipeline/test/my_project/kneaddata_result --output-prefix sample01.kneaddata --run-trim-repetitive --sequencer-source TruSeq3 -t 10 -p 10 --fastqc fastqc --remove-intermediate-output 
echo "finished kneaddata." 
/USB/meta_genomics_pipeline/bin/kraken2 --db /USB/meta_genomics_pipeline/database/kraken2_database  --threads 4  --report /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample01.kraken2.report --unclassified-out /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample01.kraken2.unclassified.fastq --output /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample01.kraken2.output /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/sample01.kneaddata.fastq 2>/USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample01.kraken2.log 
/USB/tools/Bracken-2.6.2/bracken -d /USB/meta_genomics_pipeline/database/kraken2_database -i /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample01.kraken2.report -o /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample01.S.bracken -w /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample01.S.bracken.report -r 250 -l S 
echo "finshed kraken2 and bracken on sample01" 
ktImportTaxonomy -q 2 -t 3 /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample01.kraken2.output -o /USB/meta_genomics_pipeline/test/my_project/krona_result/sample01.kraken.krona.html
/USB/meta_genomics_pipeline/bin/humann --verbose --remove-temp-output --threads 30 --input /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/sample01.kneaddata.fastq --output /USB/meta_genomics_pipeline/test/my_project/humann_result --metaphlan-options "--bowtie2db /USB/meta_genomics_pipeline/database/metaphlan3_database -t rel_ab --add_viruses "
echo "finshed humann on sample01" 
/USB/meta_genomics_pipeline/bin/humann_rename_table --input /USB/meta_genomics_pipeline/test/my_project/humann_result/sample01_genefamilies.tsv --output /USB/meta_genomics_pipeline/test/my_project/humann_result/sample01_genefamilies-names.tsv --names uniref90 
/USB/meta_genomics_pipeline/bin/humann_renorm_table --input /USB/meta_genomics_pipeline/test/my_project/humann_result/sample01_genefamilies.tsv --output /USB/meta_genomics_pipeline/test/my_project/humann_result/sample01_genefamilies-cpm.tsv --units cpm --update-snames 
/USB/meta_genomics_pipeline/bin/megahit -r /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/sample01.kneaddata.fastq -o /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01 --out-prefix sample01.megahit.final --num-cpu-threads 10 --memory 0.9 --min-contig-len 1000 
echo "finshed megahit on sample01" 
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/viral  -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.viral.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/fungi  -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.fungi.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/plasmid -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.plasmid.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/bacteria -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.bacteria.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/archaea  -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.archaea.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/UniVec_Core -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.UniVec_Core.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
cat /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.*.blastn.out > /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.cat.all.blast.out
python /USB/meta_genomics_pipeline/bin/blastout2contig_annotation.py -c /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa -b /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01/sample01.cat.all.blast.out -m /USB/meta_genomics_pipeline/database/blastn_dbs/seqid2taxid.6abfpuv.map -o /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample01.annotated.final.contigs.fa
/USB/meta_genomics_pipeline/bin/metaquast.py /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa -o /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit-quast-report 
echo "finshed metaquast on sample01" 
# 
echo "start prokka on sample01" 
prokka /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample01/sample01.megahit.final.contigs.fa --outdir /USB/meta_genomics_pipeline/test/my_project/prokka_result --force --prefix sample01 --db /home/ceid/miniconda3/envs/prokka_env/db --metagenome --kingdom Bacteria --quiet 
echo "finshed prokka on sample01" 
echo "start salmon on sample01" 
salmon index --index /USB/meta_genomics_pipeline/test/my_project/salmon_result/sample01 --transcripts /USB/meta_genomics_pipeline/test/my_project/prokka_result/sample01.ffn 
salmon quant --index /USB/meta_genomics_pipeline/test/my_project/salmon_result/sample01 --libType IU -r /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/sample01.kneaddata.fastq -o /USB/meta_genomics_pipeline/test/my_project/salmon_result/sample01 --quiet 
echo "finshed salmon on sample01" 
# 
# circos ./circos.config 
ln -sb /USB/meta_genomics_pipeline/test/test_data/HD48R4_subsample.fastq.gz /USB/meta_genomics_pipeline/test/my_project/fastqc_result/sample02.raw.fastq.gz 
echo "finshed link fastq files" 
kneaddata --input /USB/meta_genomics_pipeline/test/my_project/fastqc_result/sample02.raw.fastq.gz --reference-db /USB/meta_genomics_pipeline/database/kneaddata_database --output /USB/meta_genomics_pipeline/test/my_project/kneaddata_result --output-prefix sample02.kneaddata --run-trim-repetitive --sequencer-source TruSeq3 -t 10 -p 10 --fastqc fastqc --remove-intermediate-output 
echo "finished kneaddata." 
/USB/meta_genomics_pipeline/bin/kraken2 --db /USB/meta_genomics_pipeline/database/kraken2_database  --threads 4  --report /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample02.kraken2.report --unclassified-out /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample02.kraken2.unclassified.fastq --output /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample02.kraken2.output /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/sample02.kneaddata.fastq 2>/USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample02.kraken2.log 
/USB/tools/Bracken-2.6.2/bracken -d /USB/meta_genomics_pipeline/database/kraken2_database -i /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample02.kraken2.report -o /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample02.S.bracken -w /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample02.S.bracken.report -r 250 -l S 
echo "finshed kraken2 and bracken on sample02" 
ktImportTaxonomy -q 2 -t 3 /USB/meta_genomics_pipeline/test/my_project/kraken2_result/sample02.kraken2.output -o /USB/meta_genomics_pipeline/test/my_project/krona_result/sample02.kraken.krona.html
/USB/meta_genomics_pipeline/bin/humann --verbose --remove-temp-output --threads 30 --input /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/sample02.kneaddata.fastq --output /USB/meta_genomics_pipeline/test/my_project/humann_result --metaphlan-options "--bowtie2db /USB/meta_genomics_pipeline/database/metaphlan3_database -t rel_ab --add_viruses "
echo "finshed humann on sample02" 
/USB/meta_genomics_pipeline/bin/humann_rename_table --input /USB/meta_genomics_pipeline/test/my_project/humann_result/sample02_genefamilies.tsv --output /USB/meta_genomics_pipeline/test/my_project/humann_result/sample02_genefamilies-names.tsv --names uniref90 
/USB/meta_genomics_pipeline/bin/humann_renorm_table --input /USB/meta_genomics_pipeline/test/my_project/humann_result/sample02_genefamilies.tsv --output /USB/meta_genomics_pipeline/test/my_project/humann_result/sample02_genefamilies-cpm.tsv --units cpm --update-snames 
/USB/meta_genomics_pipeline/bin/megahit -r /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/sample02.kneaddata.fastq -o /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02 --out-prefix sample02.megahit.final --num-cpu-threads 10 --memory 0.9 --min-contig-len 1000 
echo "finshed megahit on sample02" 
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/viral  -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.viral.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/fungi  -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.fungi.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/plasmid -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.plasmid.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/bacteria -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.bacteria.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/archaea  -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.archaea.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
blastn -db /USB/meta_genomics_pipeline/database/blastn_dbs/UniVec_Core -query /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa -out /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.UniVec_Core.blastn.out -outfmt 6 -num_threads 32 -max_target_seqs 1
cat /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.*.blastn.out > /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.cat.all.blast.out
python /USB/meta_genomics_pipeline/bin/blastout2contig_annotation.py -c /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa -b /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02/sample02.cat.all.blast.out -m /USB/meta_genomics_pipeline/database/blastn_dbs/seqid2taxid.6abfpuv.map -o /USB/meta_genomics_pipeline/test/my_project/blastn_result/sample02.annotated.final.contigs.fa
/USB/meta_genomics_pipeline/bin/metaquast.py /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa -o /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit-quast-report 
echo "finshed metaquast on sample02" 
# 
echo "start prokka on sample02" 
prokka /USB/meta_genomics_pipeline/test/my_project/megahit_result/sample02/sample02.megahit.final.contigs.fa --outdir /USB/meta_genomics_pipeline/test/my_project/prokka_result --force --prefix sample02 --db /home/ceid/miniconda3/envs/prokka_env/db --metagenome --kingdom Bacteria --quiet 
echo "finshed prokka on sample02" 
echo "start salmon on sample02" 
salmon index --index /USB/meta_genomics_pipeline/test/my_project/salmon_result/sample02 --transcripts /USB/meta_genomics_pipeline/test/my_project/prokka_result/sample02.ffn 
salmon quant --index /USB/meta_genomics_pipeline/test/my_project/salmon_result/sample02 --libType IU -r /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/sample02.kneaddata.fastq -o /USB/meta_genomics_pipeline/test/my_project/salmon_result/sample02 --quiet 
echo "finshed salmon on sample02" 
# 
# circos ./circos.config 

# fastqc raw data and filgered data 
fastqc /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/*.kneaddata.fastq --outdir /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/fastqc -t 6 --quiet 
multiqc /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/fastqc --outdir /USB/meta_genomics_pipeline/test/my_project/kneaddata_result/fastqc --quiet --no-ansi 

# summary statistics: 
/USB/meta_genomics_pipeline/bin/kneaddata_read_count_table --input /USB/meta_genomics_pipeline/test/my_project/kneaddata_result --output /USB/meta_genomics_pipeline/test/my_project/summary_report/Summary.kneadata.results.xls 
/USB/meta_genomics_pipeline/bin/summary_kraken_count_table.py --input /USB/meta_genomics_pipeline/test/my_project/kraken2_result --output /USB/meta_genomics_pipeline/test/my_project/summary_report/Summary.kraken2.results.xls 
python3 /USB/meta_genomics_pipeline/bin/megahit-quest-multi-result.py --input /USB/meta_genomics_pipeline/test/my_project/megahit_result --output /USB/meta_genomics_pipeline/test/my_project/summary_report/Summary.megahit.quest.results.xls
python3 /USB/meta_genomics_pipeline/bin/prokka_multi_result.py --input /USB/meta_genomics_pipeline/test/my_project/prokka_result --output /USB/meta_genomics_pipeline/test/my_project/summary_report/Summary.Prokka.results.xls
python3 /USB/meta_genomics_pipeline/bin/gather-counts.py /USB/meta_genomics_pipeline/test/my_project/salmon_result
python3 /USB/meta_genomics_pipeline/bin/results2html.py
