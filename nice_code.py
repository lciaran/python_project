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
import flexcalc
from statistics import mean, stdev
import dictionaries


def uniprot_to_pdb(query_ID):
    """ Funtion that obtains the fasta file from Uniprot using the Uniprot ID"""
    url_link = 'https://www.uniprot.org/uniprot/' + query_ID + '.fasta'
    wget.download(url_link)

uniprot_to_pdb("P65206")

def query_info_from_fasta(fasta_file):
 """ Function that returns the query information, ID and sequence,
 from fasta file"""
 with open(fasta_file) as handle:
     for record in SeqIO.parse(handle, "fasta"):
         return (record.id, str(record.seq))

query = query_info_from_fasta("P65206.fasta")
##print (query)


def top_10_blast_idlist(fasta_file):
    """Function that performs the Blast and returns a list of the IDs from the
    top 10 results of the Blast"""
    ## performing the Blast
    cline = NcbiblastpCommandline(query=fasta_file, db="DB/PDB_db", evalue=0.00001, out= "blast_results.out", outfmt = "6 sseqid evalue")
    stdt, stdr= cline()

    ## getting the 10 best porteins IDs
    with open ("blast_results.out", "r") as file:
        list_IDs = []
        for line in file:
            ID = line[4:8]
            if len(list_IDs) < 5:
                if ID not in list_IDs:
                    list_IDs.append(ID)
    return (list_IDs)

list = top_10_blast_idlist("P65206.fasta")
##print (list)

def homologous_PDB(list_hom, query):
    """Function that obtains the PDB files of the list of top 10 homologous proteins,
    extracts the sequence and introduces them into the alignment file as well
    as the query (query has to be a tupple (id,seq))"""
    with open ("aln_input.fa", "w") as file:
        file.write(str(">" + query[0] + "\n" + query[1] +"\n"))
        for Id in list_hom:
            url_pdb = 'https://files.rcsb.org/view/' + Id + '.pdb'
            wget.download(url_pdb)
            PDB_file_path = Id + '.pdb'
            query_seqres = SeqIO.parse(PDB_file_path, 'pdb-atom')
            query_chain_id = Id.upper() + ':A'
            for chain in query_seqres:
                if chain.id == query_chain_id:
                    query_chain = chain.seq
                    file.write(str(">" + Id + "\n" + query_chain + "\n"))

homologous_PDB(list, query)

def clustalw(aln_file = "aln_input.fa"):
    """Function that performs the ClustalW alignment and converts it to fasta
    format, no input needed"""
    ## performing the clustalw
    cmd = ClustalwCommandline("clustalw", infile=aln_file, outfile="aln_output.aln")
    stdout, stderr = cmd()
    ## converting it to fasta format
    align = AlignIO.read("aln_output.aln", "clustal")
    SeqIO.write(align, "aln_output.fa", "fasta")

clustalw()

def pdb_bfactor_info(list_hom):
    """ Function that returns a dictionary of all the b-factors stored in the
    PDBs of the top proteins from the list of the blast results"""

    pdb_data = {}

    for Id in list_hom:
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
                            if line[16] != "B":
                                location += 1
                                residue = three_to_one(residue)
                                pdb_data.setdefault(Id, {}).setdefault(location, {}).setdefault(residue, bfactor)
    return(pdb_data)

dic_pdb_data= pdb_bfactor_info(list)
#print (dic_pdb_data)

def pdb_bfactor_info_normalized(pdb_data_dict):
    """Function that normalises the b-factors of the PDB and returns the
    modified dictionary"""
    for id in pdb_data_dict:
        list_bfactors = []
        for pos in pdb_data_dict[id]:
            for aa, bfactor in pdb_data_dict[id][pos].items():
                list_bfactors.append(bfactor)
        mn = mean(list_bfactors)
        std = stdev(list_bfactors)
        for pos in pdb_data_dict[id]:
            for aa, bfactor in pdb_data_dict[id][pos].items():
                pdb_data_dict[id][pos][aa] = abs((bfactor - mn) / std)
    return (pdb_data_dict)

dic_pdb_data = pdb_bfactor_info_normalized(dic_pdb_data)
##print(dic_pdb_data)

def alignment_to_dict(alignment_fasta = "aln_output.fa"):
    """Function that returns a dictionary with the postions of each aminoacid or
    gap of each protein in the fasta alignment file"""
    msa_dict = {}
    with open(alignment_fasta, "r") as fasta_msa:
        for line in SeqIO.parse(fasta_msa, "fasta"):
            position = 0
            for AA in str(line.seq):
                msa_dict.setdefault(line.id, {}).setdefault(position, AA)
                position += 1
    return (msa_dict)

dic_msa = alignment_to_dict()
##print (dic_msa)

def b_factor_dictionary(aln_dict, PDB_dict, query):
    """Function that returns the B-factor dictionary taking into account the
    alignment results and doing the mean of the values of the possible B-factors
    from the homologues"""
    ## creating the b_factor dictionary
    b_factor_dict = {}
    for id, pos in aln_dict.items():
        if id == query[0]:
            for id2, pos2 in aln_dict.items():
                aln_counter = 0
                pdb_loc = 1
                if id2 != id :
                    for aa in pos2:
                        if pos[aln_counter] == pos2[aln_counter]:
                            aminoacid = pos2[aln_counter]
                            bf_pdb = PDB_dict[id2][pdb_loc][aminoacid]
                            b_factor_dict.setdefault(aln_counter, {}).setdefault(aminoacid, []).append(bf_pdb)
                            pdb_loc = pdb_loc + 1
                        aln_counter += 1
    ## calculating the means for each position
    with open ("predicted_bfactors.txt", "w") as file:
        i = 1
        file.write(str("Position"+"\t"+"Aminoacid"+"\t"+"B-factor"+"\t"+"Type"+"\n"))
        while (i < len(query[1])):
            if i in b_factor_dict:
                for aa, list in b_factor_dict[i].items():
                    b_factor= round(sum(list) / len(list), 2)
                    if aa in dictionaries.rigid:
                        file.write(str(str(i)+"\t"+aa+"\t"+str(b_factor)+"\t"+"R"+"\n"))
                    else:
                        file.write(str(str(i)+"\t"+aa+"\t"+str(b_factor)+"\t"+"F"+"\n"))
            else:
                b_factor = flexcalc.flexcalc(query[1], i-1)
                if query[1][i-1] in dictionaries.rigid:
                    file.write(str(str(i)+"\t"+query[1][i-1]+"\t"+str(b_factor)+"\t"+"R"+"\n"))
                else:
                    file.write(str(str(i)+"\t"+query[1][i-1]+"\t"+str(b_factor)+"\t"+"F"+"\n"))
            i += 1

b_factor_dictionary(dic_msa, dic_pdb_data, query)
##print(bf_dict)
