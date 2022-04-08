import sys
import argparse
import os
from nice_code import *
from graphical_representations import *

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
				default= "./predicted_bfactors.txt",
				help="Ouput file. If not defined, the file will be named 'predicted_bfactors.txt'.")

    parser.add_argument(	'-f', '--figure',
				dest="outfigure",
				action="store",
				default= "./flexibility_plots.png",
				help="Flexibility plots file. If not defined, the file will be named 'flexibility_plots.png'.")

    parser.add_argument(	'-v', '--verbose',
				dest="verbose",
				action="store_true",
				default= False,
				help= "Print progression log to standard error.")

    options = parser.parse_args()


	# CAPTURING THE INPUT
    input_fasta = options.infile

    if len(sys.argv) == 1:
        sys.stderr.write("No input provided. Please, try again. \n")
        exit()
    else:
        if os.path.isfile(input_fasta):
            query = query_info_from_fasta(input_fasta)
        else:
            uniprot_to_pdb(input_fasta)
            input_fasta = input_fasta + ".fasta"
            query = query_info_from_fasta(input_fasta)

    if options.verbose:
	    sys.stderr.write("Query %s has been created. \n BLASTP is being carried out. \n" %query[0])

	# READING THE FASTA FILE, PERFORMING BLASTP AND GATHERING THE CLOSEST HOMOLOGOUS
    list_top10 = top_10_blast_idlist(input_fasta)

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
    outfile_file = options.outfile

    b_factor_dictionary(dic_msa, dic_pdb_data, query, outfile_file)

    if options.verbose:
	    sys.stderr.write("Program finished. Results stored in %s \n" %outfile_file)

    #OBTAINING THE RESULTS PLOT
    figure = options.outfigure

    flexibility_plots(outfile_file, figure)
