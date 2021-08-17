import sys
import os
from Bio import SeqIO


def deal_contig_file(afile):
    adict = {}
    for record in SeqIO.parse(afile,"fasta"):
        adict[str(record.id)] = str(record.seq)
    return adict

def deal_blastn_file(afile, adict,sample_name):
    name_list = []
    with open(afile,"r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            a = line.split()
            if a[0] not in name_list:
                name_list.append(a[0])
    for each in name_list:
        print(">{}_{} length={}\n{}".format(each,sample_name,len(adict[each]),adict[each]))
    
    return


blast_result_dir = sys.argv[1]
contig_result_dir = sys.argv[2]

for path,dirs,file_list in os.walk(blast_result_dir):
    for file_name in file_list:
        if file_name.endswith("blastn.corona.out"):
            file_path = os.path.join(path,file_name)
            # print(file_path)
            file_prefix = file_name.split(".")[0]
            file_prefix_path = file_prefix + '/'+ file_prefix +'.megahit.final.contigs.fa' 
            contig_name = os.path.join(contig_result_dir,file_prefix_path)
            contig_dict = deal_contig_file(contig_name)
            deal_blastn_file(file_path, contig_dict,file_prefix)
    
