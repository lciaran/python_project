import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as plt

bfactors = pd.read_csv("predicted_bfactors.txt", sep="\t", header = 0)
bfactors.head()

# sns.set_theme(style="whitegrid")

# rs = np.random.RandomState(365)
# values = rs.randn(365, 4).cumsum(axis=0)
# dates = pd.date_range("1 1 2016", periods=365, freq="D")
# data = pd.DataFrame(values, dates, columns=["A", "B", "C", "D"])
#data = data.rolling(7).mean()

#sns.lineplot(data=data, palette="tab10", linewidth=2.5)

plot = sns.scatterplot(data=bfactors, x="Position", y="B-factor")

# plot = sns.histplot(data=bfactors, x="Type", bins=30, color = "olive")

fig = plot.get_figure()
fig.savefig("output.png")