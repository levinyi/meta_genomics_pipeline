#!/usr/bin/python3
import os
import sys
import argparse
import configparser


def usage():
    print('''
    updates:
    
    20210411:   add PE data parsing.
    20210720:   optimized some parmarters.
    20200611:   created.
    '''.format())

def _argparse():
    parser = argparse.ArgumentParser(description="This is metagenomics pipeline")
    parser.add_argument('-c', '--config', action='store', dest='config', default='config.txt', help = "this is config file.")
    parser.add_argument('-l', '--list', action='store', dest='data_file', default='data.list', help = "your data list file.")
    return parser.parse_args()

def make_dir(*dir):
    for each in dir:
        if not os.path.exists(each):
            os.mkdir(each)

def read_raw_data(data_file):
    data_dict = {}
    with open(data_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            sample_name, fastq = line.split()
            #fastq1,fastq2 = fastq.split((",",";"))
            data_dict[sample_name] = os.path.abspath(fastq)
    return data_dict

def main():
    """ docstring for __main__"""
    parser = _argparse()
    cf = configparser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)
    if not cf.has_section('config'):
        sys.exit("Error: your config file is not correct.")

    # read config file:
    config_dict = {
        # script dir
        'script_dir' : os.path.abspath(cf.get('config', 'script_dir')),
        'python' : cf.get('config','python'),
        # kneaddata
        'kneaddata' : cf.get('kneaddata', 'kneaddata'),
        'kneaddata_ref': cf.get('kneaddata', 'kneaddata_ref'),
        'sequence_source': cf.get('kneaddata', 'sequence_source'),
        # kraken2
        'kraken2'       : os.path.abspath(cf.get('kraken2', 'kraken2')),
        'kraken2_db'    : os.path.abspath(cf.get('kraken2', 'kraken2_db')),
        'kraken2_threads': cf.get('kraken2','kraken2_threads'),
        'bracken'       : os.path.abspath(cf.get('braken', 'bracken')),
        'classification_level': cf.get('braken', 'classification_level'),
        'bracken_read_length' : cf.get('braken', 'bracken_read_length'),
        'kreport2mpa'   : os.path.abspath(cf.get('braken', 'kreport2mpa')),
        # humann3
        'humann' : os.path.abspath(cf.get('humann', 'humann')),
        'humann_db': os.path.abspath(cf.get('humann', 'humann_db')),
        'metaphlan3_db' :os.path.abspath(cf.get('humann', 'metaphlan3_db')),
        # blastn
        'blastn_dbs': os.path.abspath(cf.get('blastn','blastn_dbs')),
        'outfmt' : cf.get('blastn','outfmt'),
        'num_threads' : cf.get('blastn', 'num_threads'),
        # megahit
        'megahit' : os.path.abspath(cf.get('megahit', 'megahit')),
        'min-contig-len' : cf.get('megahit', 'min-contig-len'),
        # quast
         'quast' : os.path.abspath(cf.get('quast', 'quast')),
        # prokka
        'prokka' : cf.get('prokka', 'prokka'),
        # cdhit
        'cd-hit' : os.path.abspath(cf.get('cdhit','cd-hit')),
        # salmon
        'salmon' : os.path.abspath(cf.get('salmon','salmon')),
        # project name.
        'project_name' : cf.get('config','project_name'),
        }
    # print(config_dict)   # for debug.
    project_dir = os.path.abspath(".") + '/' + config_dict['project_name']
    make_dir(project_dir)
    fastqc_dir    = os.path.abspath(project_dir) + '/fastqc_result'
    kneaddata_dir = os.path.abspath(project_dir) + '/kneaddata_result'
    taxonomic_dir = os.path.abspath(project_dir) + '/kraken2_result'
    humann_dir    = os.path.abspath(project_dir) + '/humann_result'
    prokka_dir    = os.path.abspath(project_dir) + '/prokka_result'
    megahit_dir   = os.path.abspath(project_dir) + '/megahit_result'
    summary_dir   = os.path.abspath(project_dir) + '/summary_report'
    salmon_dir    = os.path.abspath(project_dir) + '/salmon_result'
    mapping_dir   = os.path.abspath(project_dir) + '/mapping_result'
    krona_dir     = os.path.abspath(project_dir) + '/krona_result'
    blastn_dir    = os.path.abspath(project_dir) + '/blastn_result'
    
    config_dict.update(
            {"project_dir" : project_dir,
              "fastqc_dir" : fastqc_dir,
           "kneaddata_dir" : kneaddata_dir,
           "taxonomic_dir" : taxonomic_dir,
               "krona_dir" : krona_dir,
              "humann_dir" : humann_dir,
             "megahit_dir" : megahit_dir, 
             "salmon_dir"  : salmon_dir,
             "prokka_dir"  : prokka_dir,
             "summary_dir" : summary_dir,
             "mapping_dir" : mapping_dir,
             "blastn_dir"  : blastn_dir,
             })
    make_dir(fastqc_dir, kneaddata_dir, taxonomic_dir, humann_dir, megahit_dir, salmon_dir, summary_dir, mapping_dir)
    make_dir(krona_dir, blastn_dir)
    print("# Create work directory")
    # print(config_dict)

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.metagenomics.sh'
    # deal fastq:
    fastq_data_dict = read_raw_data(parser.data_file)
    # print(fastq_data_dict.values()) # for debug
    ### fastq_data_dict :{'sample01': 'houmaohu_R1.fq.gz', 'sample02': '/YJ2021620_R1.fq.gz'}
    with open(shell_name, "w") as f:
        # step1: filtering.
        for each_sample in fastq_data_dict:
            print(each_sample)
            f.write("ln -sb {0} {fastqc_dir}/{1}.raw.fastq.gz \n".format(fastq_data_dict[each_sample], each_sample, **config_dict))
            f.write("echo \"finshed link fastq files\" \n")
            f.write("{kneaddata} --input {fastqc_dir}/{1}.raw.fastq.gz --reference-db {kneaddata_ref} --output {kneaddata_dir} --output-prefix {1}.kneaddata --run-trim-repetitive --sequencer-source {sequence_source} -t 10 -p 10 --fastqc fastqc --remove-intermediate-output \n".format(fastq_data_dict[each_sample], each_sample, **config_dict))
            f.write("echo \"finished kneaddata.\" \n")
        # step2 kraken2
            f.write("{kraken2} --db {kraken2_db}  --threads {kraken2_threads}  --report {taxonomic_dir}/{0}.kraken2.report --unclassified-out {taxonomic_dir}/{0}.kraken2.unclassified.fastq --output {taxonomic_dir}/{0}.kraken2.output {kneaddata_dir}/{0}.kneaddata.fastq 2>{taxonomic_dir}/{0}.kraken2.log \n".format(each_sample, **config_dict))
            f.write("{bracken} -d {kraken2_db} -i {taxonomic_dir}/{0}.kraken2.report -o {taxonomic_dir}/{0}.S.bracken -w {taxonomic_dir}/{0}.S.bracken.report -r {bracken_read_length} -l {classification_level} \n".format(each_sample, **config_dict))
            f.write("echo \"finshed kraken2 and bracken on {0}\" \n".format(each_sample))
        # krona for graph
            f.write("ktImportTaxonomy -q 2 -t 3 {taxonomic_dir}/{0}.kraken2.output -o {krona_dir}/{0}.kraken.krona.html\n".format(each_sample, **config_dict))

        # step3 humann3 : see: https://github.com/biobakery/biobakery/wiki/humann 3.1, 3.2
            f.write("{humann} --verbose --remove-temp-output --threads 30 --input {kneaddata_dir}/{0}.kneaddata.fastq --output {humann_dir} --metaphlan-options \"--bowtie2db {metaphlan3_db} -t rel_ab --add_viruses \"\n".format(each_sample, **config_dict))
            f.write("echo \"finshed humann on {0}\" \n".format(each_sample))
            f.write("{script_dir}/humann_rename_table --input {humann_dir}/{0}_genefamilies.tsv --output {humann_dir}/{0}_genefamilies-names.tsv --names uniref90 \n".format(each_sample, **config_dict))
        # Normalizing RPKs to relative abundance
            f.write("{script_dir}/humann_renorm_table --input {humann_dir}/{0}_genefamilies.tsv --output {humann_dir}/{0}_genefamilies-cpm.tsv --units cpm --update-snames \n".format(each_sample, **config_dict))
        
        # step4 assemble by megahit.
            f.write("{megahit} -r {kneaddata_dir}/{0}.kneaddata.fastq -o {megahit_dir}/{0} --out-prefix {0}.megahit.final --num-cpu-threads 10 --memory 0.9 --min-contig-len 1000 \n".format(each_sample, **config_dict))
            f.write("echo \"finshed megahit on {0}\" \n".format(each_sample))

        # blastn
        # Note: building reference index before using blastn.
            make_dir("{blastn_dir}/{0}".format(each_sample,**config_dict))
            f.write("blastn -db {blastn_dbs}/viral  -query {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -out {blastn_dir}/{0}/{0}.viral.blastn.out -outfmt {outfmt} -num_threads {num_threads} -max_target_seqs 1\n".format(each_sample, **config_dict))
            f.write("blastn -db {blastn_dbs}/fungi  -query {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -out {blastn_dir}/{0}/{0}.fungi.blastn.out -outfmt {outfmt} -num_threads {num_threads} -max_target_seqs 1\n".format(each_sample, **config_dict))
            f.write("blastn -db {blastn_dbs}/plasmid -query {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -out {blastn_dir}/{0}/{0}.plasmid.blastn.out -outfmt {outfmt} -num_threads {num_threads} -max_target_seqs 1\n".format(each_sample, **config_dict))
            f.write("blastn -db {blastn_dbs}/bacteria -query {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -out {blastn_dir}/{0}/{0}.bacteria.blastn.out -outfmt {outfmt} -num_threads {num_threads} -max_target_seqs 1\n".format(each_sample, **config_dict))
            f.write("blastn -db {blastn_dbs}/archaea  -query {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -out {blastn_dir}/{0}/{0}.archaea.blastn.out -outfmt {outfmt} -num_threads {num_threads} -max_target_seqs 1\n".format(each_sample, **config_dict))
            f.write("blastn -db {blastn_dbs}/UniVec_Core -query {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -out {blastn_dir}/{0}/{0}.UniVec_Core.blastn.out -outfmt {outfmt} -num_threads {num_threads} -max_target_seqs 1\n".format(each_sample, **config_dict))
            f.write("cat {blastn_dir}/{0}/{0}.*.blastn.out > {blastn_dir}/{0}/{0}.cat.all.blast.out\n".format(each_sample, **config_dict))
            f.write("python {script_dir}/blastout2contig_annotation.py -c {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -b {blastn_dir}/{0}/{0}.cat.all.blast.out -m {blastn_dbs}/seqid2taxid.6abfpuv.map -o {blastn_dir}/{0}.annotated.final.contigs.fa\n".format(each_sample, **config_dict))
        # mapping 
            # f.write("bowtie2-build {megahit_dir}/{0}/{0}.megahit.final.contigs.fa {megahit_dir}/{0}/{0} \n".format(each_sample, **config_dict))
            # f.write("bowtie2 --threads 10 -x {megahit_dir}/{0}/{0} -U {kneaddata_dir}/{0}.kneaddata.fastq |samtools view -@ 10 -F 4 -bS -o {mapping_dir}/{0}.bowtie2megahit.contig.bam - \n".format(each_sample, **config_dict))
            # f.write("echo \"finshed bowtie2 mapping on {0}\" \n".format(each_sample))
        # quast
            f.write("{quast} {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -o {megahit_dir}/{0}/{0}.megahit-quast-report \n".format(each_sample, **config_dict))
            f.write("echo \"finshed metaquast on {0}\" \n".format(each_sample))
            f.write("# \n".format(**config_dict))
        # prokka
            f.write("echo \"start prokka on {0}\" \n".format(each_sample))
            f.write("{prokka} {megahit_dir}/{0}/{0}.megahit.final.contigs.fa --outdir {prokka_dir} --force --prefix {0} --db /home/ceid/miniconda3/envs/prokka_env/db --metagenome --kingdom Bacteria --quiet \n".format(each_sample, **config_dict))
            f.write("echo \"finshed prokka on {0}\" \n".format(each_sample))
        # salmon
            f.write("echo \"start salmon on {0}\" \n".format(each_sample))
            f.write("salmon index --index {salmon_dir}/{0} --transcripts {prokka_dir}/{0}.ffn \n".format(each_sample, **config_dict))
            f.write("salmon quant --index {salmon_dir}/{0} --libType IU -r {kneaddata_dir}/{0}.kneaddata.fastq -o {salmon_dir}/{0} --quiet \n".format(each_sample, **config_dict))
            f.write("echo \"finshed salmon on {0}\" \n".format(each_sample))
        # circos
            f.write("# \n")
            f.write("# circos ./circos.config \n".format(**config_dict))

        # fastqc
        f.write("\n# fastqc raw data and filgered data \n")
        f.write("fastqc {kneaddata_dir}/*.kneaddata.fastq --outdir {kneaddata_dir}/fastqc -t 6 --quiet \n".format(**config_dict))
        f.write("multiqc {kneaddata_dir}/fastqc --outdir {kneaddata_dir}/fastqc --quiet --no-ansi \n".format(**config_dict))
        f.write("\n# summary statistics: \n")

        f.write("{script_dir}/kneaddata_read_count_table --input {kneaddata_dir} --output {summary_dir}/Summary.kneadata.results.xls \n".format(**config_dict))
        f.write("{script_dir}/summary_kraken_count_table.py --input {taxonomic_dir} --output {summary_dir}/Summary.kraken2.results.xls \n".format(**config_dict))
        ## f.write("{script_dir}/merge_metaphlan_tables.py {humann_dir}/*/*.tsv > ./{summary_dir}/Summary.metaphlan.results.xls \n".format(**config_dict))
        # f.write("{script_dir}/humann_join_tables --input {humann_dir} --output {summary_dir}/Summary.humann.results.xls \n".format(**config_dict))
        # Join HUMAnN3 output per sample into one table.
        # f.write("{script_dir}/humann_join_tables -s --input {humann_dir} --file_name pathabundance --output {summary_dir}/Summary.humann.pathabundance.tsv\n".format(**config_dict))
        # f.write("{script_dir}/humann_join_tables -s --input {humann_dir} --file_name pathcoverage --output {summary_dir}/Summary.humann.pathcoverage.tsv \n".format(**config_dict))
        # f.write("{script_dir}/humann_join_tables -s --input {humann_dir} --file_name genefamilies --output {summary_dir}/Summary.humann.genefamilies.tsv \n".format(**config_dict))
        # Re-normalize gene family and pathway abundances (so that all samples are in units of copies per million).
        # f.write("{script_dir}/humann_renorm_table --input {summary_dir}/Summary.humann.pathabundance.tsv --units cpm --output {summary_dir}/Summary.humann.pathabundance_cpm.tsv\n".format(**config_dict))
        # f.write("{script_dir}/humann_renorm_table --input {summary_dir}/Summary.humann.genefamilies.tsv --units cpm --output {summary_dir}/Summary.humann.genefamilies_cpm.tsv\n".format(**config_dict))
        # 
        f.write("{python} {script_dir}/megahit-quest-multi-result.py --input {megahit_dir} --output {summary_dir}/Summary.megahit.quest.results.xls\n".format(**config_dict))
        f.write("{python} {script_dir}/prokka_multi_result.py --input {prokka_dir} --output {summary_dir}/Summary.Prokka.results.xls\n".format(**config_dict))
        f.write("{python} {script_dir}/gather-counts.py {salmon_dir}\n".format(**config_dict))
        f.write("{python} {script_dir}/results2html.py\n".format(**config_dict))

    print("all finished!")

if __name__ == '__main__':
    main()
