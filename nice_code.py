import wget
from Bio import SeqIO
from sys import argv
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import numpy as np
import math
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import Polypeptide

## to obtain the fasta file for the query from the Uniprot website
query_ID = 'Q8IM15'
url_link = 'https://www.uniprot.org/uniprot/' + query_ID + '.fasta'
fasta_file = wget.download(url_link)
with open(fasta_file) as handle:
  for record in SeqIO.parse(handle, "fasta"):
    query_seq = record.seq

## performing the blastp
cline = NcbiblastpCommandline(query=fasta_file, db="DB/PDB_db", evalue=0.00001, out= "results.out", outfmt = "6 sseqid evalue")
stdt, stdr= cline()

## getting the 10 best porteins IDs
with open ("results.out", "r") as file:
    list_IDs = []
    for line in file:
        ID = line[0:4]
        if len(list_IDs) < 5:
            if ID not in list_IDs:
                list_IDs.append(ID)
#print (list_IDs)

## obtaining the PDB files and the sequences, and writing them on the alignment file
with open ("aln_input.fa", "w") as file2:
    file2.write(str(">" + query_ID + "\n" + query_seq +"\n"))
    for Id in list_IDs:
        url_pdb = 'https://files.rcsb.org/view/' + Id + '.pdb'
        pdb_file = wget.download(url_pdb)
        PDB_file_path = Id + '.pdb'
        query_seqres = SeqIO.parse(PDB_file_path, 'pdb-atom')
        query_chain_id = Id.upper() + ':A'
        for chain in query_seqres:
            if chain.id == query_chain_id:
                query_chain = chain.seq
                file2.write(str(">" + Id + "\n" + query_chain + "\n"))

## performing the clustalw
cmd = ClustalwCommandline("clustalw", infile="aln_input.fa", outfile="output.aln")
stdout, stderr = cmd()

align = AlignIO.read("output.aln", "clustal")

count = SeqIO.write(align, "output.fa", "fasta")

#b-factors from PDB:
pdb_data = {}
for Id in list_IDs:
    location = 0
    with open(Id + '.pdb') as input_pdb:
        for line in input_pdb:
            if line.startswith("ATOM"):
                chain = line[21]
                residue = line[17:20]
                atom = line[13:16]
                bfactor = float(line[61:66])
                if 'A' in chain:
                    if 'CA' in atom:
                        location += 1
                        residue = three_to_one(residue)
                        pdb_data.setdefault(Id, {}).setdefault(residue, {}).setdefault(location, bfactor)

#first 3AA
# trip_dict = {}
# for id in pdb_data:
#   triplet = ''
#   for residue in pdb_data[id]:
#     for location in pdb_data[id][residue]:
#       triplet += list(pdb_data[id].keys())[list(pdb_data[id][residue].keys()).index(location)]
#   trip_dict[id] = [triplet[0:3]]
#   trip_dict[id].append(triplet[-3:])


#from fasta to b-factors
msa_dict = {}
with open("output.fa", "r") as fasta_msa:
  for line in SeqIO.parse(fasta_msa, "fasta"):
    msa_dict[line.id] = str(line.seq)

# b_factor = {}
# for ID1, fa_seq1 in msa_dict.items():
#     for ID2, fa_seq2 in msa_dict.items():
#         if ID1 != ID2:
#             t1 = triplet[0]
#             t2 = triplet[1]
#             i1 = fa_seq.index(t1)
#             i2 = fa_seq.index(t2) + 2
#             msa_reduced_dict[ID1] = fa_seq[i1:i2]

bf_dict = {}

seq1 = msa_dict[query_ID]
position = 0

while position < len(seq1):
    for id2, seq2 in msa_dict.items():
        if query_ID != id2:
            for AA1 in seq1:
                for AA2 in seq2:
                    if AA1 == AA2 and seq1[position] == seq2[position]:
                        bf_list = []
                        location = 0 + len(bf_list)
                        bf_pdb = pdb_data[id2][AA2][location]
                        bf_list.append(bf_pdb)
                        print(bf_list)
                    else:
                        pass
                    #else:
                        #article function (bf_calculated)
                        #bf_list.append(bf_calculated)
#                 bf_query = sum(bf_list) / len(bf_list)
#                 bf_dict[AA1] = bf_query
            position = position + 1

# print (bf_dict)
