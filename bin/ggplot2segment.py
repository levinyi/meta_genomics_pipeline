import sys
from Bio import SeqIO


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


reference = sys.argv[1]
blast_result = sys.argv[2]

ref_y = 10
print("refname\tx1\tx2\ty1\ty2")
for record in SeqIO.parse(reference,"fasta"):
    print("ref\t1\t{}\t{}\t{}".format(len(record.seq),ref_y,ref_y))

adict = {}
with open(blast_result,"r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        line = line.rstrip("\n")
        a = line.split()
        addtwodimdict(adict,a[0],a[6],a[7])

##########just for debug
# import json
# print(json.dumps(adict,indent=4))
#####################################

for k,v in adict.items():
    ref_y = ref_y-1
    for x,y in v.items():
        ref_y = ref_y -0.5
        print("{}\t{}\t{}\t{}\t{}".format(k,x,y,ref_y,ref_y))


