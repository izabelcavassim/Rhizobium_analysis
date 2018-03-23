import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

# initiliaze a dataframe with index and column names
idf = pd.DataFrame.from_items([('A', [1, 2, 3]), ('B', [4, 5, 6]), ('C', [10, 20, 30]), ('D', [14, 15, 16])], orient='index', columns=['x', 'y', 'z'])

# Plot the clustermap which will be a figure by itself
cax = sns.clustermap(idf, col_cluster=False, row_cluster=True)

# Get the column dendrogram axis
cax_col_dend_ax = cax.ax_col_dendrogram.axes

# Plot the boxplot on the column dendrogram axis
# I still need to figure out how to show the axis for this boxplot
idf.plot(kind='box', ax=cax_col_dend_ax)

# Show the plot
plt.show()