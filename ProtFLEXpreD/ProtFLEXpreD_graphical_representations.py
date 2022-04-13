import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as plt

def flexibility_plots(predicted_bfactors_file, figure_file, bioP_bfactors_file = "./Intermediary/flexibility_bioP.txt"):
    '''Function that transforms two text files with Position and B-factors as columns into a dataframe to plot them.
    It obtains three plots that relate B-factors and their positions. Two of them represent the results obtained by the program.
    The other one represents the results obtained using the biopython module.'''
    bfactors = pd.read_csv(predicted_bfactors_file, sep="\t", header = 0)
    bfactors_bioP = pd.read_csv(bioP_bfactors_file, sep="\t", header = 0)
    fig, axs = plt.pyplot.subplots(3, 1, figsize=(50, 50))
    axs[0].set_title('Flexibility plot', fontsize=60)
    axs[1].set_title('Flexibility plot calculated with biopython', fontsize=60)
    axs[2].set_title('B-factors coloured by property', fontsize=60)
    sns.lineplot(data=bfactors, x="Position", y="B-factor", color="olive", linewidth=3.5, ax=axs[0])
    sns.lineplot(data=bfactors_bioP, x="Position", y="B-factor", color="purple", linewidth=3.5, ax=axs[1])
    sns.scatterplot(data=bfactors, x="Position", y="B-factor", hue="Type", s=300, ax=axs[2])
    axs[2].legend(loc=2, markerscale=5., scatterpoints=1, fontsize = 50)
    axs[0].set_xlabel(xlabel = '', fontsize=55)
    axs[0].set_ylabel(ylabel = '', fontsize=55)
    axs[1].set_ylabel(ylabel = 'B-factors', fontsize=55)
    axs[1].set_xlabel(xlabel = '', fontsize=55)
    axs[2].set_ylabel(ylabel = '', fontsize=55)
    axs[2].set_xlabel(xlabel = 'Position', fontsize=55)
    axs[0].tick_params(axis='both', which='major', labelsize=50)
    axs[1].tick_params(axis='both', which='major', labelsize=50)
    axs[2].tick_params(axis='both', which='major', labelsize=50)
    fig = fig.get_figure()
    fig.savefig(figure_file)

flexibility_plots('./Results/P65206.txt', './Results/P65206.png')
