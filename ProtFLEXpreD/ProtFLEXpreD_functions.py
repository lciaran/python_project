import sys
import re
import wget
from Bio import SeqIO
from Bio import AlignIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Polypeptide import three_to_one
from statistics import mean, stdev
import dictionaries


def uniprot_to_pdb(query_ID):
    """ Funtion that obtains the fasta file from Uniprot using the Uniprot ID"""
    try:
        url_link = 'https://www.uniprot.org/uniprot/' + query_ID + '.fasta'
        wget.download(url_link, './Downloads')
    except:
        sys.stderr.write("Please enter a valid Uniprot ID.")
        exit()

def query_info_from_fasta(fasta_file):
    """ Function that returns the query information, ID and sequence,
    from fasta file"""
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return (record.id, str(record.seq))

def top_10_blast_idlist(fasta_file):
    """Function that performs the Blast and returns a list of the IDs from the
    top 10 results of the Blast"""
    ## performing the Blast
    cline = NcbiblastpCommandline(query=fasta_file, db="DB_uniprot/uniprot_sprot.fasta", evalue=0.00001, out= "./Alignment/blast_results.out", outfmt = "6 sseqid evalue")
    stdt, stdr= cline()

    ## getting the 10 best porteins IDs
    with open ("./Alignment/blast_results.out", "r") as file:
        list_IDs = []
        for line in file:
            ID = line[3:9]
            if len(list_IDs) < 10:
                if ID not in list_IDs:
                    list_IDs.append(ID)
    return (list_IDs)

def homologous_PDB(list_hom, query):
    """Function that obtains the PDB files of the list of top 10 homologous proteins,
    extracts the sequence and introduces them into the alignment file as well
    as the query (query has to be a tupple (id,seq))"""
    pdb_data = {}
    with open ("./Alignment/aln_input.fa", "w") as file:
        try:
            file.write(str(">" + query[0] + "\n" + query[1] +"\n"))
        except TypeError:
            sys.stderr.write("Please enter a file in fasta format.\n")
            exit()
        for Id in list_hom:
            try:
                url_pdb = "https://alphafold.ebi.ac.uk/files/AF-" + Id + "-F1-model_v2.pdb"
                wget.download(url_pdb, './Downloads')
            except:
                continue
            PDB_file_path = "./Downloads/AF-" + Id + "-F1-model_v2.pdb"
            sequence = ''
            try:
                with open (PDB_file_path, "r") as pdb_file:
                    for line in pdb_file:
                        if line.startswith("ATOM"):
                            chain = line[21]
                            residue = line[17:20]
                            atom = line[13:16]
                            location = int(line[22:27])
                            bfactor = float(line[61:66])
                            if 'A' in chain:
                                if 'CA' in atom:
                                    residue = three_to_one(residue)
                                    sequence = sequence + residue
                                    pdb_data.setdefault(Id, {}).setdefault(location, {}).setdefault(residue, bfactor)
                file.write(str(">" + Id + "\n" + sequence + "\n"))
            except ValueError:
                pass
    return pdb_data

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
                try:
                    pdb_data_dict[id][pos][aa] = abs((bfactor - mn) / std)
                except ZeroDivisionError:
                    std = 1
                    pdb_data_dict[id][pos][aa] = abs((bfactor - mn) / std)
    return (pdb_data_dict)

def clustalw(aln_file = "./Alignment/aln_input.fa"):
    """Function that performs the ClustalW alignment and converts it to fasta
    format, no input needed"""
    ## performing the clustalw
    cmd = ClustalwCommandline("clustalw", infile=aln_file, outfile="./Alignment/aln_output.aln")
    stdout, stderr = cmd()
    ## converting it to fasta format
    align = AlignIO.read("./Alignment/aln_output.aln", "clustal")
    SeqIO.write(align, "./Alignment/aln_output.fa", "fasta")

def alignment_to_dict(alignment_fasta = "./Alignment/aln_output.fa"):
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

def flexcalc(protein, r):
    '''Function that calculates the residues b-factors according to their neighbors'''             
    b_query =  None
    query = protein[r]
    if r == 0:
        b_query = dictionaries.b_factor[query]
    elif r == len(protein) - 1:
        b_query = dictionaries.b_factor[query]
    else:
        if protein[r - 1] in dictionaries.rigid:
            if protein[r + 1] in dictionaries.rigid:
                b_query = dictionaries.b_factor_rr[query]
            elif protein[r + 1] in dictionaries.flex:
                b_query = dictionaries.b_factor_rf[query]
        elif protein[r - 1] in dictionaries.flex:
            if protein[r + 1] in dictionaries.rigid:
                b_query = dictionaries.b_factor_rf[query]
            elif protein[r + 1] in dictionaries.flex:
                b_query = dictionaries.b_factor_ff[query]
    return b_query

def b_factor_dictionary(aln_dict, PDB_dict, query, output_file):
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
                            if pos[aln_counter] != "-":
                                aminoacid = pos2[aln_counter]
                                bf_pdb = PDB_dict[id2][pdb_loc][aminoacid]
                                b_factor_dict.setdefault(aln_counter, {}).setdefault(aminoacid, []).append(bf_pdb)
                                pdb_loc = pdb_loc + 1
                        elif re.search("[A-Z]", pos2[aln_counter]):
                            pdb_loc = pdb_loc + 1
                        aln_counter += 1
    ## calculating the means for each position
    with open (output_file, "w") as file:
        i = 0
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
                b_factor = flexcalc(query[1], i)
                if query[1][i] in dictionaries.rigid:
                    file.write(str(str(i)+"\t"+query[1][i]+"\t"+str(b_factor)+"\t"+"R"+"\n"))
                else:
                    file.write(str(str(i)+"\t"+query[1][i]+"\t"+str(b_factor)+"\t"+"F"+"\n"))
            i += 1

def flex_bioP(sequence, flexibility_file = "./Alignment/flexibility_bioP.txt"):
    '''Function that returns a file with the b-factor of each position. It uses the flexibility
    function provided by biopython library.'''
    analysed_seq = ProteinAnalysis(sequence)
    flexibility = analysed_seq.flexibility()
    with open(flexibility_file, "w") as flexi_file:
        flexi_file.write("Position" + "\t" + "B-factor" + "\n")
        p = 0
        for i in flexibility: 
            flexi_file.write(str(p) + "\t" + str(round(i, 4)) + "\n")
            p += 1
