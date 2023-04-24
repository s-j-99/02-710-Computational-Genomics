import networkx as nx
import pandas as pd
from grandiso import find_motifs
import matplotlib.pyplot as plt
import random
import numpy as np


# Load the data and view the columns
all_data = pd.read_csv(
    "./glioblastoma/BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS-4.4.220.tab3.csv")
print(all_data.columns)

# Filter out only human interactions
human_ppi = all_data[(all_data['Organism Name Interactor A'] == "Homo sapiens") & (
    all_data['Organism Name Interactor B'] == "Homo sapiens")]

print(
    f"glioblastoma human interactions has {human_ppi.shape[0]} rows")

# Create Undirected Graph
Graphtype = nx.Graph()
g_glio_human = nx.from_pandas_edgelist(
    human_ppi, create_using=Graphtype, source="Official Symbol Interactor A", target='Official Symbol Interactor B')
print(g_glio_human)
# Look at frequencies of node degrees
deg_hist = nx.degree_histogram(g_glio_human)
plt.bar(x=[i for i in range(20)], height=deg_hist[:20])
plt.title("Degree frequencies in glioblastoma network")
# plt.show()
# plt.clf()

# Find average node degree in network
av_deg = 0
counts = 0
for idx, i in enumerate(deg_hist):
    counts += i
    av_deg += (idx + 1) * i
av_deg = av_deg / counts

print(
    f"glioblastoma human interactions has an average of {av_deg} edges per node")

density = nx.density(g_glio_human)

print(
    f"glioblastoma human interactions has a density of {density}")

clust = nx.clustering(g_glio_human)
av_clust_coef = np.array(list(clust.values())).mean()
print(
    f"glioblastoma human interactions has an average clustering coefficient of {av_clust_coef}")

# Find the number of connected components in entire network
n_conn_comp = nx.number_connected_components(g_glio_human)
print(
    f"glioblastoma human interactions has {n_conn_comp} connected components")
