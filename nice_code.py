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
query_ID = 'P06401'
url_link = 'https://www.uniprot.org/uniprot/' + query_ID + '.fasta'
fasta_file = wget.download(url_link)
with open(fasta_file) as handle:
  for record in SeqIO.parse(handle, "fasta"):
    query_seq= record.seq

# ## performing the blastp
cline = NcbiblastpCommandline(query=fasta_file, db="DB/PDB_db", evalue=0.00001, out= "results.out", outfmt = "6 sseqid evalue")
stdt, stdr= cline()

## getting the 10 best porteins IDs
with open ("results.out", "r") as file:
    list_IDs = []
    for line in file:
        ID = line[4:8]
        if len(list_IDs) < 2:
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
                        pdb_data.setdefault(Id, {}).setdefault(location, {}).setdefault(residue, bfactor)

##print (pdb_data)

#from fasta to b-factors
msa_dict = {}
with open("output.fa", "r") as fasta_msa:
  for line in SeqIO.parse(fasta_msa, "fasta"):
      position = 0
      for AA in str(line.seq):
        msa_dict.setdefault(line.id, {}).setdefault(position, AA)
        position += 1

## print(msa_dict)

## creating the b_factor dictionary
b_factor_dict = {}

for id, seq in msa_dict.items():
    if id == query_ID:
        bf_list = []
        for id2, seq2 in msa_dict.items():
            aln_counter = 0
            pdb_loc = 1
            if id2 != id :
                for aa in seq:
                    if seq[aln_counter] == seq2[aln_counter]:
                        aa = seq2[aln_counter]
                        bf_pdb = pdb_data[id2][pdb_loc][aa]
                        b_factor_dict.setdefault(aln_counter, {}).setdefault(aa, []).append(bf_pdb)
                        pdb_loc = pdb_loc + 1
                    aln_counter += 1

## print(b_factor_dict)
i = 1
while (i < len(query_seq)):
    if i in b_factor_dict.keys():
        for aa, list in b_factor_dict[i].items():
            b_factor_dict[i][aa] = round(sum(list) / len(list), 2)
    #else:
        #setdefault()
    i += 1

## print(b_factor_dict)
