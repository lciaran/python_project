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
    query_seq= record.seq

## performing the blastp
cline = NcbiblastpCommandline(query=fasta_file, db="DB/PDB_db", evalue=0.00001, out= "results.out", outfmt = "6 sseqid evalue")
stdt, stdr= cline()

## getting the 10 best porteins IDs
with open ("results.out", "r") as file:
    list_IDs = []
    for line in file:
        ID = line[4:8]
        if len(list_IDs) < 10:
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
        query_seqres = SeqIO.parse(PDB_file_path, 'pdb-seqres')
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
    with open(Id + '.pdb') as input_pdb:
        for line in input_pdb:
            if line.startswith("ATOM"):
                chain = line[21]
                residue = line[17:20]
                location = line[23:26]
                atom = line[13:16]
                bfactor = line[61:66]
                if 'A' in chain:
                    if 'CA' in atom:
                        residue = three_to_one(residue)
                        pdb_data.setdefault(Id, {}).setdefault(residue, {}).setdefault(location, []).append(bfactor)

##print (pdb_data)
