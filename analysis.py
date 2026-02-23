# %% [markdown]
# #libraries used

# %%
import pandas as pd
import glob
import os

# %% [markdown]
# #Finding all quant.sf files

# %%
files = glob.glob("*/quant.sf")
print(files)
count_dict = {}

# %% [markdown]
# #Making the count_matrix table

# %%
for f in files:
    sample = os.path.basename(os.path.dirname(f))
    df = pd.read_csv(f, sep="\t")
    count_dict[sample] = df["NumReads"].values


import numpy as np

count_df = pd.DataFrame(count_dict)
gene_ids = pd.read_csv("ALLgene_ids.csv", sep="\t")["Name"]


count_df.insert(0, "Gene_ID", gene_ids)
count_df.columns = ["Gene_ID", "sample1", "sample2", "sample3", "sample4"]   

count_df = count_df.set_index("Gene_ID")
count_df = count_df.groupby("Gene_ID").sum()

numeric_cols = ["sample1", "sample2", "sample3", "sample4"]
count_df = count_df[count_df[numeric_cols].sum(axis = 1) > 0]

count_df = count_df.T
count_df




# %% [markdown]
# #Building the metadata

# %%
metadata_dict = {
    "sample" : ["sample1", "sample2", "sample3", "sample4"],
    "condition" : ["control", "control", "treated", "treated"]
}


metadata = pd.DataFrame(metadata_dict)
metadata

# %%
#Saving csv files

# %%
count_df.to_csv("counts_matrix.csv", index=False)
metadata.to_csv("metadata.csv",index=False)

# %% [markdown]
# #Rebuilding the table with gene_ids

# %% [markdown]
# #Merging 2 files

# %%
count_df = pd.read_csv("Counts_matrix2.csv", sep= ",")
count_df = count_df.reset_index().rename(columns={"Name" : "Gene_ID"})
if "index" in count_df.columns:
    count_df = count_df.drop(columns=["index"])

count_df

# %%
count_df.to_csv("Final_Counts_matrix.csv", index=False)

# %% [markdown]
# #Aggregation

# %%
count_df = pd.read_csv('Final_Counts_matrix.csv', sep= ",")
count_df = count_df.groupby("Gene_ID").sum().reset_index()
count_df

# %%
count_df = count_df.set_index("Gene_ID")
count_df

# %%
count_df = count_df[count_df[numeric_cols].sum(axis = 1) > 0]
numeric_cols = ["sample1", "sample2", "sample3", "sample4"]
count_df

# %%
count_df = count_df.T
count_df

# %%
metadata = metadata.set_index("sample")
metadata

# %% [markdown]
# #DE analysis

# %% [markdown]
# #Converting table values to int
# #Creating the dds object
# #deseq method

# %%
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

count_cols = count_df.columns
count_df[count_cols] = count_df[count_cols].astype(int)

dds = DeseqDataSet(counts = count_df, 
metadata = metadata,
design = "condition")


dds.deseq2()
type(dds)



# %%
stat_res = DeseqStats(dds, n_cpus=8, contrast=("condition", "control", "treated"))
stat_res.summary()



# %%
res = stat_res.results_df
res = res.reset_index()
res = res.rename(columns={"index" : "Gene_ID"})
res


# %% [markdown]
# #mapping gene 

# %%
gtf_map = pd.read_csv("gtf.mapping.csv", sep= " ", names=["Gene_ID", "Gene_Symbol"])
res = res.drop(columns=[c for c in res.columns if "Gene_Symbol" in c], errors="ignore")


res = res.reset_index()
res = res.merge(gtf_map, on="Gene_ID", how="left")
res= res.set_index("Gene_ID")
res


# %%

volcano = res.copy()
volcano["sig"] = "NS"   

volcano["neg_log10_padj"] = -np.log10(volcano["padj"])
volcano.loc[(volcano["padj"] < 0.05) & (volcano["log2FoldChange"] > 1), "sig"] = "Up"
volcano.loc[(volcano["padj"] < 0.05) & (volcano["log2FoldChange"] < -1), "sig"] = "Down"

colors = {"NS": "grey", "Up": "red", "Down": "blue"}

plt.scatter(
    volcano["log2FoldChange"],
    volcano["neg_log10_padj"],
    c=volcano["sig"].map(colors),
    s=5,
    alpha=0.95
)

plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
plt.axvline(1, color="black", linestyle="--", linewidth=1)
plt.axvline(-1, color="black", linestyle="--", linewidth=1)

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 Adjusted p-value")
plt.title("Volcano plot: Budesonide-treated vs untreated asthma cells")

plt.show()



import numpy as np
import matplotlib.pyplot as plt

#top 30 most diffrentially expressed genes
up = res[(res["log2FoldChange"]) > 2 & (res["padj"] < 0.05)]
down = res[(res["log2FoldChange"]) < -2 & (res["padj"] < 0.05)]

top30_up = up.sort_values("padj").head(30)
top30_down = down.sort_values("padj").head(30)
print(top30_up)
top30_down





# %%
top30_up


# %% [markdown]
# #heatmap plot

# %%
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler


dds.layers["log1p"] = np.log1p(dds.layers["normed_counts"])
dds_down = dds[:, top30_down.index]
grapher = pd.DataFrame(dds_down.layers["log1p"].T, index=dds_down.var_names, columns=dds_down.obs_names)
id_to_symbol = top30_down["Gene_Symbol"]   
#print(id_to_symbol) 
grapher.index = grapher.index.map(id_to_symbol)

grapher = grapher.dropna(axis=0)                 # remove genes with no symbol
grapher = grapher.groupby(grapher.index).mean()  # merge duplicate symbols safely
grapher = grapher.loc[grapher.std(axis=1) > 0]

print(grapher)

sns.clustermap(grapher, z_score=0, cmap= "RdYlBu_r")





# %%
dds.layers["log1p"] = np.log1p(dds.layers["normed_counts"])
dds_up = dds[:, top30_up.index]
grapher = pd.DataFrame(dds_up.layers["log1p"].T, index=dds_up.var_names, columns=dds_down.obs_names)
id_to_symbol = top30_up["Gene_Symbol"]   
#print(id_to_symbol) 
grapher.index = grapher.index.map(id_to_symbol)

grapher = grapher.dropna(axis=0)                 # remove genes with no symbol
grapher = grapher.groupby(grapher.index).mean()  # merge duplicate symbols safely
grapher = grapher.loc[grapher.std(axis=1) > 0]

print(grapher)

g = sns.clustermap(grapher, z_score=0, cmap= "RdYlBu_r")

g.ax_heatmap.set_title(
    "Heatmap of downregulated genes\nBudesonide-treated vs untreated asthma cells",
    fontsize=14,
    pad=30
)

g.ax_heatmap.set_ylabel("Gene symbols")
g



# %%
abs_sig = res[(abs(res["log2FoldChange"]) > 3) & (res["padj"] < 0.5)]
#up = res[(res["log2FoldChange"]) > 2 & (res["padj"] < 0.05)]
abs_sig = abs_sig.sort_values("padj").head(30)


dds.layers["log1p"] = np.log1p(dds.layers["normed_counts"])
dds_abs = dds[:, abs_sig.index]
grapher = pd.DataFrame(dds_abs.layers["log1p"].T, index=dds_abs.var_names, columns=dds_abs.obs_names)
id_to_symbol = abs_sig["Gene_Symbol"]   
#print(id_to_symbol) 
grapher.index = grapher.index.map(id_to_symbol)

grapher = grapher.dropna(axis=0)                 # remove genes with no symbol
grapher = grapher.groupby(grapher.index).mean()  # merge duplicate symbols safely
grapher = grapher.loc[grapher.std(axis=1) > 0]


print(abs_sig)

g = sns.clustermap(grapher, z_score=0, cmap= "RdYlBu_r")
g.ax_heatmap.set_title(
    "Heatmap of most diffrentially expressed genes\nBudesonide-treated vs untreated asthma cells",
    fontsize=14,
    pad=30
)

g.ax_heatmap.set_ylabel("Gene symbols")
g



