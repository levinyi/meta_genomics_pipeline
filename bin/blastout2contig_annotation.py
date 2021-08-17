#!/home/ceid/miniconda3/bin/python3
import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
import argparse
import decimal


def usage():
	print("""
		usage:
			python {} -c <contig_file> -b <blastn_file> -m <id2map_file> -o <output_file>

	20210728: create.

	these files are used for test:
	"/USB/meta_genomics_pipeline/database/blastn_dbs/seqid2taxid.6abfpuv.map"
	"/USB/metagenomics/analysis/my_project/megahit_result/houmaohu/houmaohu.megahit.final.contigs.fa"
	"/USB/metagenomics/analysis/my_project/blastn_result/houmaohu/houmaohu.cat.all.blastn.out"
		""".format(sys.argv[0]))


def _argparse():
    parser = argparse.ArgumentParser(description="This is xxx")
    parser.add_argument('-c', '--contig_file', action='store', dest='contig_file', help = "contig file from Megahit result.")
    parser.add_argument('-b', '--blastn_file', action='store', dest='blastn_file', help = "blastn result file.")
    parser.add_argument('-m', '--id2map_file', action='store', dest='id2map_file', help = "id to map file.")
    parser.add_argument('-o', '--output_file', action='store', dest='output_file', help = "output file,fasta format.")
    return parser.parse_args()


def main():
    """ docstring for __main__"""
    parser = _argparse()

    ######################
    # deal with id2map file return a dict.

    taxonomy_dict = {}
    seq2len_gc_dict = {}
    with open(parser.id2map_file, "r") as f:
        for line in f:
            line = line.rstrip("\n").lstrip(">")
            a, length, gc = line.split("\t")
            # print(a)
            b = a.split()
            name = b[0]
            # print(name)
            desc = " ".join(b[1:])
            # print(desc)
            taxonomy_dict[name] = desc # taxonomy_dict = {"kraken:taxid|60847|NZ_CP048739.1" : "Mycobacteroides xxx chromosome, complete genome"}
            seq2len_gc_dict[name] = length +' '+ gc

    ################
    # deal blast file return a dataframe.
    blast_tab = pd.read_table(parser.blastn_file, sep="\t", header=None)
    new_blast_tab = blast_tab.drop_duplicates([0])
    new_blast_tab.columns=["a", "b", "c","d","e","f","g","h","i","j","k","l"]
    new_blast_tab = new_blast_tab.set_index("a")

    ######################
    # deal with contig fasta file, output a new fasta file with detailed descriptions.
    output_table = open(str(parser.output_file)[:-3]+'.table.txt', "w")
    output_table.write("contig_name\tannotation\tref_name\tref_gc\tcontig_len\tcontig_gc\tGenome_coverage\n")
    output_handle = open(parser.output_file, "w")
    for record in SeqIO.parse(parser.contig_file,"fasta"):
        if record.id in new_blast_tab.index:
            annotation1 = new_blast_tab.loc[record.id, "b"]  # record.id=k141_208, annotation1 = kraken:taxid|60847|NZ_CP048739.1
            annotation2 = taxonomy_dict[annotation1]  # annotation2="Mycobacteroides xxx chromosome, complete genome"
            contig_leng = len(record.seq)
            contig_GC = GC(record.seq)
            record.description = record.description + " contig_GC="+"%.2f"%contig_GC+" ref"+seq2len_gc_dict[annotation1] +" "+ annotation1 + " " + annotation2 
            SeqIO.write(record, output_handle, "fasta")
            ref_len,ref_gc = seq2len_gc_dict[annotation1].split()
            output_table.write("{}\t{}\t{}\t{}\t{}\t{:.2f}\t{:.3%}\n".format(record.id, annotation2, ref_len, ref_gc, contig_leng, contig_GC,float(contig_leng)/int(ref_len)))
        else:
            record.description = record.description + ' None of Blastn result'
            # print("None blastn result")
    output_handle.close()
    output_table.close()


if __name__ == '__main__':
	main()
