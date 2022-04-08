import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt

def flexibility_plots(predicted_bfactors_file, figure_file):
    bfactors = pd.read_csv(predicted_bfactors_file, sep="\t", header = 0)
    bfactors.head()
    fig, axs = plt.pyplot.subplots(2, 1, figsize=(25, 10))
    axs[0].set_title('Flexibility plot', fontsize=15)
    axs[1].set_title('B-factors coloured by property', fontsize=15)
    sns.lineplot(data=bfactors, x="Position", y="B-factor", color="olive", linewidth=2.5, ax=axs[0])
    sns.scatterplot(data=bfactors, x="Position", y="B-factor", hue="Type", ax=axs[1])
    fig = fig.get_figure()
    fig.savefig(figure_file)