# 02-710-Computational-Genomics
Network Motif Finding

This project aims to implement some motif-finding algorithms on a biological dataset. The dataset used here is the Glioblastoma Protein-Protein Interactions data compiled by BioGRID.

Two types of motif-search strategies are employed:
1. Using predefined motifs: Some known biological motifs are Feed-Forward Loops (FFLs), 3-Cycles, BiFans, Single-Input Module (SIM) and Multiple Input Module (MIM). In particular, we search for FFLs, 3-Cycles and BiFans in a directed version of the graph.
2. Searching for motifs de novo: When interesting motifs are not known, we can use existing motif structures in the network to query the whole network. However, this can be computationally expensive and is generally considered NP-complete. The pattern-growing and pruning method is used here to grow size-k motifs from scratch and query them in the networks.

Motifs are counted using the subgraph census method rather than an exhaustive search. This involves taking a large number of subsamples of the original network and searching for motifs within those subgraphs to create a frequency distribution of motif counts. This should take on a Gaussian form which can be compared with random networks using a standard Z-test. We use randomly generated Barabasi-Albert (BA) and Erdos-Renyi networks here to create "background" frequency distributions.

Finally, biological significance of the highest occurring size-k motifs in the undirected network is assessed. This is done by choosing the motifs with the highest frequency, taking all the proteins found in every occurrence, and passing this list to the Gene Ontology (GO) Resource. A list of all proteins found in the glioblastoma network is also used as a reference list to calculate the enrichment of terms in specific motifs. We find highly enriched and statistically significant Biological Processes (BP) relating to cell signaling, protein modifications and metabolism, as well as PANTHER pathways relating to epithelial-mesenchymal transition (EMT) and brain-specific processes such as axon guidance and dopamine signaling. 

The code used in this project, the results and the written report are available in this repository.
