 if __name__ == "__main__":


	parser = argparse.ArgumentParser(description="Given DNA FASTA file(s), calculate the sequence length and molecular weight of their corresponding translated proteins. Outputs the output sorted by sequence length")

	parser.add_argument( 	'-i', '--input',
				dest= "infile",
				action= "store",
				default= "./",
				help="Input FASTA file or directory containing fasta file")

	parser.add_argument(	'-o', '--output',
				dest="outfile",
				action="store",
				default= None,
				help="Ouput file. If not defined, it prints output to standard output.")

	parser.add_argument(	'-v', '--verbose',
				dest="verbose",
				action="store_true",
				default= False,
				help= "Print progression log to standard error")

	parser.add_argument( 	'-p', '--pattern',
				dest="pattern",
				action="store",
				default=None,
				help="Regular expression pattern to search in the translated sequences")

	parser.add_argument(	'-r', '--random',
        	                dest="random_output_num",
				action="store",
				default= None,
				type= int,
				help="Random how many sequence to print")

	options = parser.parse_args()


	# CAPTURING THE INPUT FILE(s)
	input_path = options.infile

	if os.path.isdir(input_path):
		list_files = [ os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith(".fa") or f.endswith(".fasta") ]
	elif os.path.isfile(input_path):
		list_files = [ input_path ]
	else:
		list_files = []

	if options.verbose:
		sys.stderr.write("%d FASTA files found.\n" %len(list_files))


	# READING THE SEQUENCE FILES AND STORE THEM IN A LIST
	output_list = []

	if options.pattern is not None:
		pattern_re = re.compile(options.pattern)


	for filename in list_files:
		for dna in FASTA_iterator(filename, DNASequence):
			prot = dna.translate()
			if options.pattern is None:
				output_list.append( ( prot.get_identifier(), len(prot), prot.get_mw()) )
			elif pattern_re.search(prot.get_sequence()):		
				output_list.append( ( prot.get_identifier(), len(prot), prot.get_mw()) )
			
		if options.verbose:
			sys.stderr.write("%s finished.\n" %filename)

	if options.verbose:
		sys.stderr.write("%s sequences found.\n" %len(output_list))

	if options.outfile is None:
		out_fd = sys.stdout
	else:
		if options.outfile.endswith(".gz"):
			out_fd = gzip.open(options.outfile, "wt")
		else:
			out_fd = open(options.outfile, "w")

	# WRITE OUTPUT
	if options.random_output_num is not None:
		if options.random_output_num > len(output_list):
			sys.stderr.write("Number of required sequences to be printed is greater than the number of available sequences\n")
			sys.exit(1)

		output_to_print = random.sample(output_list, options.random_output_num)
		
	else:
		output_to_print = output_list

	# SORT THE SEQUENCES BY THE LENGTH
	if options.verbose:
		sys.stderr.write("Sorting the sequences...\n")
	output_to_print.sort(key=lambda x: x[1], reverse=True)
	if options.verbose:
		sys.stderr.write("Sort process finished.\n")

	for identifier, length, mw in output_to_print:
		out_fd.write("%s\t%s\t%s\n" %(	identifier,
						length,
						mw))

	out_fd.close()

	if options.verbose:
		sys.stderr.write("Program finished correctly.")
