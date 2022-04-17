import sys
import os
import argparse
from ProtFLEXpreD.PFD_functions import *
from ProtFLEXpreD.PFD_representations import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Given the query protein FASTA file or Uniprot ID, calculate the protein sequence flexibility (B-factors). The output is a file containing the sequence aminoacids B-factors and two plots representing the protein's flexibility")

    parser.add_argument( 	'-i', '--input',
				dest= "infile",
				action= "store",
				default= "./",
				help="Input protein FASTA file or ID from Uniprot")

    parser.add_argument(	'-o', '--output',
				dest="outfile",
				action="store",
				default= "./Results/predicted_bfactors.txt",
				help="Ouput file. If not defined, the file will be named 'predicted_bfactors.txt'.")

    parser.add_argument(	'-f', '--figure',
				dest="outfigure",
				action="store",
				default= "./Results/flexibility_plots.png",
				help="Flexibility plots file. If not defined, the file will be named 'flexibility_plots.png'.")

    parser.add_argument(	'-v', '--verbose',
				dest="verbose",
				action="store_true",
				default= False,
				help= "Print progression log to standard error.")

    options = parser.parse_args()

    # CREATING DIRECTORIES
    os.makedirs('./Results', exist_ok=True)

    #CREATING THE PDB DATABASE IF DOES NOT EXIST
    db_files =["./DB_pdb/PDB.phr", "./DB_pdb/PDB.pin", "./DB_pdb/PDB.psq", "./DB_pdb/pdb_seqres.txt"]
    for A in db_files:
        if os.path.exists(A):
            continue
        else:
            sys.stderr.write("Downloading the PDB database.\n")
            database()

	# CAPTURING THE INPUT, THE OUTPUT AND THE FIGURE
    input_fasta = options.infile
    outfile_file = options.outfile
    figure = options.outfigure

    if len(sys.argv) == 1:
        sys.stderr.write("No input provided. Please, try again. \n")
        exit()
    else:
        if os.path.isfile(input_fasta):
            query = query_info_from_fasta(input_fasta)
        else:
            uniprot_to_pdb(input_fasta)
            input_fasta = "./Downloads/" + input_fasta + ".fasta"
            query = query_info_from_fasta(input_fasta)

    try:
        if options.verbose:
            sys.stderr.write("Query %s has been created. \n BLASTP is being carried out. \n" %query[0])
    except TypeError:
        sys.stderr.write("Please enter a file in fasta format.\n")
        exit()

	# READING THE FASTA FILE, PERFORMING BLASTP AND GATHERING THE CLOSEST HOMOLOGOUS
    list_top10 = top_10_blast_idlist(input_fasta)

    if len(list_top10) == 0:
        #ALTERNATIVE TO A PROTEIN WITHOUT HOMOLOGOUS
        if options.verbose:
            sys.stderr.write("%d homologous proteins have been found. \n Getting B-factors according to their position. \n " %len(list_top10))
        p = 0
        with open (outfile_file, "w") as file:
            file.write(str("Position"+"\t"+"Aminoacid"+"\t"+"B-factor"+"\t"+"Type"+"\n"))
            for p in query[1]:
                b_factor = flexcalc(query[1], p)
                if query[1][p] in dictionaries.rigid:
                    file.write(str(str(p)+"\t"+query[1][p]+"\t"+str(b_factor)+"\t"+"R"+"\n"))
                else:
                    file.write(str(str(p)+"\t"+query[1][p]+"\t"+str(b_factor)+"\t"+"F"+"\n"))
                p += 1
        flex_bioP(query[1])
        flexibility_plots(outfile_file, figure)
        if options.verbose:
            sys.stderr.write("Program finished. Results stored in %s \n" %outfile_file)
        exit()

    if options.verbose:
	    sys.stderr.write("%d homologous proteins have been found. \n Getting B-factors from the homologous proteins PDB files. \n" %len(list_top10))
    
    #CREATING THE DICTIONARY WITH THE INFORMATION OF THE HOMOLOGOUS PROTEINS PDBs.
    dic_pdb_data = homologous_PDB(list_top10, query)
    dic_pdb_data = pdb_bfactor_info_normalized(dic_pdb_data)

    if options.verbose:
	    sys.stderr.write("Information gathered. \n Performing clustalW and storing results. \n")

    #PERFORMING MSA AND SAVING ALIGNMENT INFORMATION IN A DICTIONARY
    clustalw()
    dic_msa = alignment_to_dict()
    
    if options.verbose:
        sys.stderr.write("ClustalW finished. \n Obtaining B-factors of each aminoacid in query sequence. \n")

    #OBTAINING B-FACTORS OF EACH AMINOACID OF THE QUERY SEQUENCE AND STORING IT INTO THE OUTPUT FILE

    b_factor_dictionary(dic_msa, dic_pdb_data, query, outfile_file)

    flex_bioP(query[1])

    if options.verbose:
	    sys.stderr.write("Program finished. Results stored in %s \n" %outfile_file)

    #OBTAINING THE RESULTS PLOT

    flexibility_plots(outfile_file, figure)
