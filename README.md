# 02-710-Computational-Genomics
Network Motif Finding

This project aims to implement some motif-finding algorithms on a biological dataset. The dataset used here is the Glioblastoma Protein-Protein Interactions data compiled by BioGRID.
Two types of motif-search strategies are employed:
1. Using predefined motifs: Some known biological motifs are Feed-Forward Loops (FFLs), BiFans, Single-Input Module (SIM) and Multiple Input Module (MIM). These motifs are searched in a directed version of the graph.
2. Searching for motifs de novo: When interesting motifs are not known, we can use existing motif structures in the network to query the whole network. However, this can be computationally expensive and is generally considered NP-complete.
