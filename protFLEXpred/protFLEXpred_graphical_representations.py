import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt

def flexibility_plots(predicted_bfactors_file, figure_file, bioP_bfactors_file = "./Alignment/flexibility_bioP.txt"):
    '''Function that transforms two text files with Position and B-factors as columns into a dataframe to plot them.
    It obtains three plots that relate B-factors and their positions. Two of them represent the results obtained by the program.
    The other one represents the results obtained using the biopython module.'''
    bfactors = pd.read_csv(predicted_bfactors_file, sep="\t", header = 0)
    bfactors_bioP = pd.read_csv(bioP_bfactors_file, sep="\t", header = 0)
    fig, axs = plt.pyplot.subplots(3, 1, figsize=(40, 20))
    axs[0].set_title('Flexibility plot', fontsize=15)
    axs[1].set_title('B-factors coloured by property', fontsize=15)
    axs[2].set_title('Flexibility plot calculated with biopython', fontsize=15)
    sns.lineplot(data=bfactors, x="Position", y="B-factor", color="olive", linewidth=2.5, ax=axs[0])
    sns.scatterplot(data=bfactors, x="Position", y="B-factor", hue="Type", ax=axs[1])
    sns.lineplot(data=bfactors_bioP, x="Position", y="B-factor", color="purple", linewidth=2.5, ax=axs[2])
    fig = fig.get_figure()
    fig.savefig(figure_file)
