import sys
import argparse
import pandas as pd
from Bio import SeqIO

def _argparse():
    parser = argparse.ArgumentParser(description="This is xxx")
    parser.add_argument('-f', '--fasta_file', action='store', dest='fasta_file', help = "fasta file")
    parser.add_argument('-k', '--key_words',   action='store', dest='key_words', help = "key words")
    parser.add_argument('-o', '--output_file', action='store', dest='output_file', help = "output file fasta format.")
    return parser.parse_args()

def main():
    parser = _argparse()
    key_words = parser.key_words
    # print(key_words)
    key_words_seq_count = 0
    key_words_seq_id = []

    fasta_file = parser.fasta_file
    sample_name = fasta_file.split(".")[0]
    with open(parser.output_file,"w") as f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if key_words in record.description:
                key_words_seq_count += 1
                key_words_seq_id.append(record.id)
                record.id = sample_name + "_" + record.id
                SeqIO.write(record, f, "fasta")
    print("{}\t{}\t{}".format(sample_name,key_words_seq_count, ",".join(key_words_seq_id)))
    # print("done")
if __name__ == '__main__':
    main()