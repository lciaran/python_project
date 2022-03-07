from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio import SeqIO
# output = commands.getoutput("clustalw --version")
#         if "not found" not in output and "CLUSTAL" in output and "Multiple Sequence Alignments" in output:
#             clustalw_exe = "clustalw"
#
# if not clustalw_exe:
#     raise MissingExternalDependencyError("Install clustalw or clustalw2 if you want to use it from Biopython.")

cmd = ClustalwCommandline("clustalw", infile="/home/laura/PYT/Exercises/2.fa", outfile="test.aln")

stdout, stderr = cmd()

align = AlignIO.read("test.aln", "clustal")

print(align)

count = SeqIO.write(align, "example.faa", "fasta")
