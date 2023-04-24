import networkx as nx
import pandas as pd
from grandiso import find_motifs
import matplotlib.pyplot as plt
import random


# Load the data and view the columns
all_data = pd.read_csv(
    "./glioblastoma/BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS-4.4.220.tab3.csv")
print(all_data.columns)

# Filter out only human interactions
human_ppi = all_data[(all_data['Organism Name Interactor A'] == "Homo sapiens") & (
    all_data['Organism Name Interactor B'] == "Homo sapiens")]

print(f"glioblastoma human interactions has {human_ppi.shape[0]} rows")

# Create Undirected Graph
Graphtype = nx.Graph()
g_glio_human = nx.from_pandas_edgelist(
    human_ppi, create_using=Graphtype, source="Official Symbol Interactor A", target='Official Symbol Interactor B')


def enumerate_subgraphs(graph, simulation_number=500, subgraph_size=100, motif_size=5):
    '''
    Generate all possible subgraphs of size=motif_size, and search for repeats in rest of graph. Avoid checking the same motif twice.
    '''
    # Simulate 500 samplings of a subgraph of size 100 from the original network
    simulation_dict = {}  # Store result of each simulation
    for i in range(simulation_number):
        # Create a random subgraph of the network
        random_subgraph_edges = random.sample(
            list(g_glio_human.edges()), subgraph_size)
        subgraph = nx.Graph()
        subgraph.add_edges_from(random_subgraph_edges)
        # Store motifs and edges already found
        motif_list = []
        edges = []
        motif_count_dict = {}  # store the motif and its count
        num_subgraph_edges = len(subgraph.edges())
        index = 0
        # terminate loop when all edges in the subgraph have been used in motifs
        while len(edges) != num_subgraph_edges:

            if index // 5 == 0:
                print("Iteration number: ", index)

            # select random size=motif_size motif
            random_motif_node_start = random.sample(
                list(subgraph.nodes()), 1)

            motif = CreateMotif(
                subgraph, random_motif_node_start, motif_size, motif_list)

            # Check if this motif was used before
            is_motif_new = CheckMotifRepeat(motif, motif_list)

            if is_motif_new and nx.is_connected(motif):

                # If new, keep in motif list
                motif_list.append(motif)
                # Add the unique edges from the motif to edge list
                edges = AppendEdges(edges, motif.edges())
                # Count the number of this motif
                matches = find_motifs(motif, subgraph)
                # Store the result with the index of the search
                motif_count_dict[index] = [motif, len(matches), matches]
                index += 1
            else:
                continue
        simulation_dict[i] = motif_count_dict
    return simulation_dict


# def CreateMotif(subgraph, random_motif_node_start, motif_size, motif_list):
#     '''
#     Creates a motif from a given node by accessing its neighbours in the subgraph. Checks if the motif is repeated.
#     '''
#     if len(random_motif_node_start) == motif_size:
#         # if we found a motif of correct size, check if the motif was already created before
#         if CheckMotifRepeat(motif, motif_list):
#             return motif
#         else:
#             # if found before, create a new motif from scratch
#             random_motif_node_start = random.sample(
#                 list(subgraph.nodes()), 1)
#             return CreateMotif(subgraph, random_motif_node_start, motif_size, motif_list, motif=nx.Graph())
#     else:
#         nbr_list = {}
#         # find all neighbours of nodes already in the motif
#         for start_node in random_motif_node_start.nodes():
#             neighbours = nx.neighbors(subgraph, start_node)
#             nbr_list[start_node] = neighbours
#         edge_list = motif.edges()
#         pick_random_edge = random.sample(neighbours)
#     return


def CheckMotifRepeat(query, motif_list):
    '''
    This checks if a query motif is found within a list of motifs. Returns true if the motif is NOT present.
    '''
    result = True
    for motif in motif_list:
        if nx.is_isomorphic(motif, query):
            result = False
            break
    return result


def AppendEdges(query, edge_list):
    '''
    This appends an edge from query into edge_list, only if it wasn't already there.
    '''
    num_edges = len(query)
    for edge in query:
        print(edge)
        if edge not in edge_list:
            edge_list.append(edge)
    return edge_list


print(enumerate_subgraphs(g_glio_human, simulation_number=1))
