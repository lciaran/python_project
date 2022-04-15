**ProtFLEXpreD**
=================

*Laura Ciaran Alfano, Neus Pou Amengual*

*MSc in Bioinformatics for Health Sciences*

*Project of SBI - PYT*

*2021 - 2022*

## **Introduction**

ProtFLEXpreD is a python package to predict *b-factors* of each protein aminoacid, which are values used to predict the protein flexibility.

The package works as follows. Firstly, **BlastP** compares the protein query to the *PDB* database to gather its top 10 homologous proteins. Secondly, *b-factors* of the homologues are obtained and normalized from the **PDB** PDB files. Thirdly, *Multiple Sequence Alignment* is performed with **ClustalW**. Fourthly, *b-factors* from homologous proteins are assigned to the aligned aminoacids by computing the mean between them, and the *b-factors* obtained, according to their neighbours, [by computation](https://onlinelibrary.wiley.com/doi/full/10.1110/ps.0236203) are assigned to the non-aligned aminoacids. Finally, a text file containing the *b-factor* of each aminoacid of the query sequence and three plots are created to represent the flexibility results.

## **Initializating the program**

Execute `setup.py` with the following commands to install all required packages. It may be possible that you have to execute it as root.

```{.sh}
python3 setup.py sdist
python3 setup.py build
python3 setup.py install
```
Then, make sure that you have the following scripts `dictionaries.py`, `__init__.py`, `ProtFLEXpreD.py`, `ProtFLEXpreD_functions.py`, `ProtFLEXpreD_graphical_representations.py` and the database folder `DB_pdb` in your directory.

## **Running the program**

Execute `ProtFLEXpreD.py` with the following command to run the program.

```{.sh}
python3 ProtFLEXpreD.py -i XXX
```

Where you have to write your **input** `-i` argument replacing `XXX` by the **Uniprot ID** or the **fasta** file of your ptotein query. You can also add the following arguments:

- `-o` : output file. It is saved as default in `/Results/predicted_bfactors.txt`.

- `-f`: output figure. It is saved as default in `/Results/flexibility_plots.png`.

- `-v`: verbose. It is `False` as default.

If you want to know more, write the following command-line.

```{.sh}
python3 ProtFLEXpreD.py -h
```

## **ProtFLEXpreD Results**

After running the program successfully, you will find three new folders in your directory:

* `Downloads` : where you can find all downloaded files from **PDB** or **Uniprot**.

* `Intermediary` : where you can find **BlastP** results, **ClustalW** alignment files and the predicted flexibility calculated by the [*biopython* package](https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html).

* `Results` : where you can find a file with the predicted *b-factors* of your query and representation plots of your results.

## **Examples**

Below you will find some examples to better understand to how execute the program.

* Executing the program while printing the progression log. The input is an uniprot ID. Output and output figure arguments are not indicated.

```{.sh}
python3 ProtFLEXpreD.py -i P65206 -v
```

* Executing the program without printing the progression log. The input is a fasta file and the output is a text file. Figure argument is not indicated.

```{.sh}
python3 ProtFLEXpreD.py -i P65206.fasta -o ./Results/P65206.txt
```

* Executing the program while printing the progression log. All possible arguments are indicated.

```{.sh}
python3 ProtFLEXpreD.py -i P65206 -o ./Results/P65206.txt -f ./Results/P65206.png -v
```

If you have any doubt, do not hesitate to contact us, laura.ciaran01@estudiant.upf.edu and neus.pou01@estudiant.upf.edu.
