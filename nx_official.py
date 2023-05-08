import random
import networkx as nx
import pandas as pd
import numpy as np
from grandiso import find_motifs
import matplotlib.pyplot as plt
import scipy.stats

# Load data and filter by experiment type (just as example)

all_data = pd.read_csv(
    "glioblastoma/BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS-4.4.220.tab3.csv")

all_data_h = all_data[(all_data['Organism Name Interactor A'] == "Homo sapiens") & (
    all_data['Organism Name Interactor B'] == "Homo sapiens")]

# Example code to predefine variouos motifs in directed graphs

# Feed-Forward Loop (FFL)

ffl = nx.DiGraph([(1, 2), (2, 3), (1, 3)])

# 3-Cycle

tcc = nx.DiGraph([(1, 2), (2, 3), (3, 1)])

# Single Input Module (SIM)


def make_sim(root_node, list_targets):
    '''
    nodes can be passed as both int or str types.
    root_node: single identifier for root node.
    list_targets: an array of identifiers for targets for the root node.
    '''
    # create graph object
    g = nx.DiGraph()
    # create root node
    g.add_node(root_node)
    # create target nodes
    targets = [list_targets[i] for i in range(len(list_targets))]
    g.add_nodes_from(targets)
    # create root-to-target edges
    edges = [(root_node, list_targets[i]) for i in range(len(list_targets))]
    g.add_edges_from(edges)
    return g

# Multiple Input Module (MIM)


def make_mim(list_roots, list_targets_per_root):
    '''
    nodes can be passed as both int or str types.
    list_roots: an array identifiers for each root node.
    list_targets_per_root: an array of arrays of identifiers for targets for each root node. targets can be shared by multiple root nodes.
    '''
    if len(list_roots) != len(list_targets_per_root):
        raise ValueError(
            "Number of targets given does not match number of roots.")
    else:
        g = nx.DiGraph()
        for i in range(len(list_roots)):
            # create a SIM for each root based on the list of targets given
            sim_make = make_sim(list_roots[i], list_targets_per_root[i])
            g.add_nodes_from(sim_make)
            g.add_edges_from(sim_make.edges)
    return g

# BiFan


def make_bifan(list_roots, list_targets_per_root):
    '''
    nodes can be passed as both int or str types.
    list_roots: an array identifiers for each root node. Only two are allowed for bifan.
    list_targets_per_root: an array of arrays of identifiers for targets for each root node. Only two are allowed for bifan.
    '''
    if (len(list_roots) != 2) or (len(list_targets_per_root) != 2):
        raise ValueError(
            "Input to bifan should have exactly 2 root nodes and 2 target nodes.")
    else:
        g = make_mim(
            list_roots, [list_targets_per_root, list_targets_per_root])
        return g


sim_test = make_sim("A", ["B", "C"])

mim_test = make_mim(["A", "B", "C"], [["D", "E", "F"], [
    "E", "F", "G"], ["E", "G", "H", "I"]])

bifan_test = make_bifan(["A", "B"], ["C", "D"])


'''
#### Uncomment below to visualize the motifs
Drawing code
nx.draw(sim_test, with_labels=True, font_weight='bold')
plt.show()
plt.clf()

nx.draw(mim_test, with_labels=True, font_weight='bold')
plt.show()
plt.clf()

nx.draw(bifan_test, with_labels=True, font_weight='bold')
plt.show()
plt.clf()

'''


def calc_average_degree(g):
    degrees = nx.degree(g)
    total = 0
    for node, degree in degrees:
        total += degree
    return total / len(degrees)


'''
Uncomment below for optimizing random graph parameters.
'''

# N = 250
# k_values = [100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300]
# iterator = 0
# ba_random_values = np.zeros(shape=(110,))
# er_random_values = np.zeros(shape=(110,))
# for k in range(30, 102, 10):
#     for m in range(1, 92, 1):
#         if k <= m:
#             continue
#         p = m / 1000.0 + 0.1
#         this_ba_graph = nx.barabasi_albert_graph(k, m).to_directed()
#         this_er_graph = nx.erdos_renyi_graph(k, p).to_directed()
#         ba_motifs = find_motifs(bifan_test, this_ba_graph)
#         ba_average_degree = calc_average_degree(this_ba_graph)
#         er_average_degree = calc_average_degree(this_ba_graph)
#         print("k: " + str(k))
#         print("m: " + str(m))
#         print(len(ba_motifs))
#         print("avg: " + str(ba_average_degree))
#         ba_random_values[iterator] = len(ba_motifs)
#         er_motifs = find_motifs(bifan_test, this_er_graph)
#         # print(len(er_motifs))
#         #print("avg: " + str(er_average_degree))
#         er_random_values[iterator] = len(er_motifs)
#         iterator += 1

# plt.scatter(k_values, ba_random_values[:12])
# plt.show()

# Load the glioblastoma data into a directed network

Graphtype = nx.DiGraph()
g_glio_all = nx.from_pandas_edgelist(all_data, source="Official Symbol Interactor A",
                                     target='Official Symbol Interactor B', create_using=Graphtype)
g_glio_human = nx.from_pandas_edgelist(all_data_h, source="Official Symbol Interactor A",
                                       target='Official Symbol Interactor B', create_using=Graphtype)
print("All: ", len(all_data.index))
print("Human: ", len(all_data_h.index))
print("All data: ", g_glio_all)
print("Human data: ", g_glio_human)

'''
This part finds the number of BiFan motifs in the 3 types of graphs.
'''
print("Running simulations for BiFan Motifs")
k = 30
m = 2
p = m / 1000.0 + 0.1
num_simulations = 100


this_ba_graph = nx.barabasi_albert_graph(k, m).to_directed()
this_er_graph = nx.erdos_renyi_graph(k, p).to_directed()

nx.draw(this_ba_graph, pos=nx.spring_layout(this_ba_graph))
plt.show()
plt.clf()
nx.draw(this_er_graph, pos=nx.spring_layout(this_er_graph))
plt.show()
plt.clf()

random_ba_motifs = np.zeros(shape=(num_simulations,))
random_er_motifs = np.zeros(shape=(num_simulations,))
glioblastoma_motifs = np.zeros(shape=(num_simulations,))

for simulation in range(num_simulations):

    this_ba_graph = nx.barabasi_albert_graph(k, m).to_directed()

    ba_motifs = find_motifs(bifan_test, this_ba_graph)
    random_ba_motifs[simulation] = len(ba_motifs)

    p = m / 1000.0 + 0.1
    this_er_graph = nx.erdos_renyi_graph(k, p).to_directed()

    er_motifs = find_motifs(bifan_test, this_er_graph)
    random_er_motifs[simulation] = len(er_motifs)

    random_sample_edges = random.sample(list(g_glio_human.edges), 3000)
    G_sample = nx.DiGraph()
    G_sample.add_edges_from(random_sample_edges)

    these_glioblastoma_motifs = find_motifs(bifan_test, G_sample)
    glioblastoma_motifs[simulation] = len(these_glioblastoma_motifs)

# approximately normalize
glioblastoma_motifs /= 10

print("Plotting graphs for BiFan motifs")

fig, axes = plt.subplots(3, 1, sharex=True)

plt.suptitle("Frequency of BiFan Motifs in Each Graph Type")

axes[0].hist(random_ba_motifs)
axes[0].set_title(
    "Number of BiFan Motifs in Random Barbasi Albert Graph", fontsize=10)

axes[1].hist(random_er_motifs)
axes[1].set_title(
    "Number of BiFan Motifs in Random Erdos Renyi Graph", fontsize=10)

axes[2].hist(glioblastoma_motifs)
axes[2].set_title(
    "Number of BiFan Motifs Across Random Subsamples of Glioblastoma Graph", fontsize=10)
plt.xlim(0, max(max(glioblastoma_motifs), max(
    random_er_motifs), max(random_ba_motifs)))
plt.tight_layout()
plt.show()

print("Running statistical tests.")
mean_glioblastoma = np.mean(glioblastoma_motifs)
mean_ba = np.mean(random_ba_motifs)
mean_er = np.mean(random_er_motifs)

std_glioblastoma = np.std(glioblastoma_motifs)
std_ba = np.std(random_ba_motifs)
std_er = np.std(random_er_motifs)

Z_ba = (mean_glioblastoma - mean_ba) / std_ba
Z_er = (mean_glioblastoma - mean_er) / std_er

p_value_ba = scipy.stats.norm.sf(abs(Z_ba))

p_value_er = scipy.stats.norm.sf(abs(Z_er))

print("### Results for BiFan Motif ###")
print("MEANS: Barabasi-Albert: {ba} || Erdos-Renyi: {er} || Glioblastoma: {gb}".format(
    ba=mean_ba, er=mean_er, gb=mean_glioblastoma))
print("SD: Barabasi-Albert: {ba} || Erdos-Renyi: {er} || Glioblastoma: {gb}".format(
    ba=std_ba, er=std_er, gb=std_glioblastoma))
print(
    "Z-SCORE: Glioblastoma vs Barabasi-Albert: {ba} || Glioblastoma vs Erdos-Renyi: {er}".format(ba=Z_ba, er=Z_er))
print("P-VALUE: Glioblastoma vs Barabasi-Albert: {ba} || Glioblastoma vs Erdos-Renyi: {er}".format(
    ba=p_value_ba, er=p_value_er))

'''
This part finds the number of FFL motifs in the 3 types of graphs.
'''

print("Running simulations for FFL Motifs")

k = 30
m = 2
num_simulations = 100
random_ba_motifs = np.zeros(shape=(num_simulations,))
random_er_motifs = np.zeros(shape=(num_simulations,))
glioblastoma_motifs = np.zeros(shape=(num_simulations,))
for simulation in range(num_simulations):

    this_ba_graph = nx.barabasi_albert_graph(k, m).to_directed()

    ba_motifs = find_motifs(ffl, this_ba_graph)
    random_ba_motifs[simulation] = len(ba_motifs)

    p = m / 1000.0 + 0.1
    this_er_graph = nx.erdos_renyi_graph(k, p).to_directed()

    er_motifs = find_motifs(ffl, this_er_graph)
    random_er_motifs[simulation] = len(er_motifs)

    random_sample_edges = random.sample(list(g_glio_human.edges), 3000)
    G_sample = nx.DiGraph()
    G_sample.add_edges_from(random_sample_edges)

    these_glioblastoma_motifs = find_motifs(ffl, G_sample)
    glioblastoma_motifs[simulation] = len(these_glioblastoma_motifs)


# approximately normalize
glioblastoma_motifs /= 10

print("Plotting graphs for FFL motifs")

fig, axes = plt.subplots(3, 1, sharex=True)

plt.suptitle("Frequency of FFL Motifs in Each Graph Type")

axes[0].hist(random_ba_motifs)
axes[0].set_title(
    "Number of FFL Motifs in Random Barbasi Albert Graph", fontsize=10)

axes[1].hist(random_er_motifs)
axes[1].set_title(
    "Number of FFL Motifs in Random Erdos Renyi Graph", fontsize=10)

axes[2].hist(glioblastoma_motifs)
axes[2].set_title(
    "Number of FFL Motifs Across Random Subsamples of Glioblastoma Graph", fontsize=10)
plt.xlim(0, 250)
plt.tight_layout()
plt.show()

print("Running statistical tests.")
mean_glioblastoma = np.mean(glioblastoma_motifs)
mean_ba = np.mean(random_ba_motifs)
mean_er = np.mean(random_er_motifs)
std_glioblastoma = np.std(glioblastoma_motifs)
std_ba = np.std(random_ba_motifs)
std_er = np.std(random_er_motifs)

Z_ba = (mean_glioblastoma - mean_ba) / std_ba
Z_er = (mean_glioblastoma - mean_er) / std_er

p_value_ba = scipy.stats.norm.sf(abs(Z_ba))

p_value_er = scipy.stats.norm.sf(abs(Z_er))

print("### Results for FFL Motif ###")
print("MEANS: Barabasi-Albert: {ba} || Erdos-Renyi: {er} || Glioblastoma: {gb}".format(
    ba=mean_ba, er=mean_er, gb=mean_glioblastoma))
print("SD: Barabasi-Albert: {ba} || Erdos-Renyi: {er} || Glioblastoma: {gb}".format(
    ba=std_ba, er=std_er, gb=std_glioblastoma))
print(
    "Z-SCORE: Glioblastoma vs Barabasi-Albert: {ba} || Glioblastoma vs Erdos-Renyi: {er}".format(ba=Z_ba, er=Z_er))
print("P-VALUE: Glioblastoma vs Barabasi-Albert: {ba} || Glioblastoma vs Erdos-Renyi: {er}".format(
    ba=p_value_ba, er=p_value_er))


'''
This part finds the number of 3-cycle motifs in the 3 types of graphs.
'''

print("Running simulations for 3-Cycle Motifs")

k = 30
m = 2
num_simulations = 100
random_ba_motifs = np.zeros(shape=(num_simulations,))
random_er_motifs = np.zeros(shape=(num_simulations,))
glioblastoma_motifs = np.zeros(shape=(num_simulations,))
for simulation in range(num_simulations):

    this_ba_graph = nx.barabasi_albert_graph(k, m).to_directed()

    ba_motifs = find_motifs(tcc, this_ba_graph)
    random_ba_motifs[simulation] = len(ba_motifs)

    p = m / 1000.0 + 0.1
    this_er_graph = nx.erdos_renyi_graph(k, p).to_directed()

    er_motifs = find_motifs(tcc, this_er_graph)
    random_er_motifs[simulation] = len(er_motifs)

    random_sample_edges = random.sample(list(g_glio_human.edges), 3000)
    G_sample = nx.DiGraph()
    G_sample.add_edges_from(random_sample_edges)

    these_glioblastoma_motifs = find_motifs(tcc, G_sample)
    glioblastoma_motifs[simulation] = len(these_glioblastoma_motifs)


# approximately normalize
glioblastoma_motifs /= 10

print("Plotting graphs for 3-Cycle motifs")

fig, axes = plt.subplots(3, 1, sharex=True)

plt.suptitle("Frequency of 3-Cycle Motifs in Each Graph Type")

axes[0].hist(random_ba_motifs)
axes[0].set_title(
    "Number of 3-Cycle Motifs in Random Barbasi Albert Graph", fontsize=10)

axes[1].hist(random_er_motifs)
axes[1].set_title(
    "Number of 3-Cycle Motifs in Random Erdos Renyi Graph", fontsize=10)

axes[2].hist(glioblastoma_motifs)
axes[2].set_title(
    "Number of 3-Cycle Motifs Across Random Subsamples of Glioblastoma Graph", fontsize=10)
plt.xlim(0, 250)
plt.tight_layout()
plt.show()

print("Running statistical tests.")
mean_glioblastoma = np.mean(glioblastoma_motifs)
mean_ba = np.mean(random_ba_motifs)
mean_er = np.mean(random_er_motifs)
std_glioblastoma = np.std(glioblastoma_motifs)
std_ba = np.std(random_ba_motifs)
std_er = np.std(random_er_motifs)

Z_ba = (mean_glioblastoma - mean_ba) / std_ba
Z_er = (mean_glioblastoma - mean_er) / std_er

p_value_ba = scipy.stats.norm.sf(abs(Z_ba))

p_value_er = scipy.stats.norm.sf(abs(Z_er))

print("### Results for 3-Cycle Motif ###")
print("MEANS: Barabasi-Albert: {ba} || Erdos-Renyi: {er} || Glioblastoma: {gb}".format(
    ba=mean_ba, er=mean_er, gb=mean_glioblastoma))
print("SD: Barabasi-Albert: {ba} || Erdos-Renyi: {er} || Glioblastoma: {gb}".format(
    ba=std_ba, er=std_er, gb=std_glioblastoma))
print(
    "Z-SCORE: Glioblastoma vs Barabasi-Albert: {ba} || Glioblastoma vs Erdos-Renyi: {er}".format(ba=Z_ba, er=Z_er))
print("P-VALUE: Glioblastoma vs Barabasi-Albert: {ba} || Glioblastoma vs Erdos-Renyi: {er}".format(
    ba=p_value_ba, er=p_value_er))
