**ProtFLEXpreD**
=================

*Laura Ciaran Alfano, Neus Pou Amengual*

*MSc in Bioinformatics for Health Sciences*

*Project of SBI - PYT*

*2021 - 2022*

## **Introduction**

ProtFLEXpreD is a python package to predict *b-factors* of each protein aminoacid, which are values used to predict the protein flexibility. 

The package works as follows. Firstly, **BlastP** compares the protein query to the *Uniprot* database to gather its top 10 homologous proteins. Secondly, *b-factors* of the homologues are obtained and normalized from the **Alphafold** PDB files. Thirdly, *Multiple Sequence Alignment* is performed with **ClustalW**. Fourthly, *b-factors* from homologous proteins are assigned to the aligned aminoacids by computing the mean between them, and the *b-factors* obtained [according to their neighbours by computation](https://www.polarmicrobes.org/protein-flexibility-calculation-with-python/) are assigned to the non-aligned aminoacids. Finally, a text file containing the *b-factor* of each aminoacid of the query sequence and three plots are created to represent the flexibility results.

## **Initializating program**

Execute `setup.py` by the following command to install all required packages. It is possible that you have to execute it as root.

```{.sh}
python3 setup.py install
```
Then, make sure that you have the following scripts `dictionaries.py`, `__init__.py`, `ProtFLEXpreD.py`, `ProtFLEXpreD_functions.py`, `ProtFLEXpreD_graphical_representations.py` and the database folder `DB_Uniprot` in your directory.

## **Running the program**

Execute `ProtFLEXpreD.py` by the following command to run the program.

```{.sh}
python3 ProtFLEXpreD.py -i XXX
```

Where you have to write your **input** `-i` argument replacing `XXX` by the **Uniprot ID** or the **fasta** file of your target. You can also add the following arguments:

- `-o` : output file. It is saved as default in `/Results/predicted_bfactors.txt`.

- `-f`: output figure. It is saved as default in `/Results/flexibility_plots.png`.

- `-v`: verbose. It is `False` as default.

If you want to know more, write the following command-line.

```{.sh}
python3 ProtFLEXpreD.py -h
```

## **ProtFLEXpreD Results**

After running the program successfully, you will find three new folders in your directory:

* `Downloads` : where you can find all downloaded files from **Alphafold** or **Uniprot**.

* `Intermediary` : where you can find **BlastP** results, **ClustalW** alignment files and the predicted flexibility calculated by the *biopython* package.

* `Results` : where you can find a file with the predicted *b-factors* of your query and representative plots of your results.

If you have any doubt, do not hesitate to contact with us, laura.ciaran01@estudiant.upf.edu and neus.pou01@estudiant.upf.edu.
