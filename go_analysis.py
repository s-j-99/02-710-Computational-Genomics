import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# Load PPI network
all_data = pd.read_csv(
    "glioblastoma/BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS-4.4.220.tab3.csv")

# Filter out only human interactors
all_data_h = all_data[(all_data['Organism Name Interactor A'] == "Homo sapiens") & (
    all_data['Organism Name Interactor B'] == "Homo sapiens")]

# Create a non-repeating list of proteins in the entire network
int_a = list(all_data_h["Official Symbol Interactor A"].values)
int_b = list(all_data_h["Official Symbol Interactor B"].values)
int_a.extend(int_b)
all_int = list(set(int_a))

# Save the list of proteins to a file
with open("all_proteins.txt", "w") as file:
    for i in all_int:
        file.writelines(i + "\n")

# Concatenate files containing proteins in every motif occurrence for size-4 motif 1

path_1 = "./size4_shape1"

# Concatenate every row
unique1 = []
with open("size4_shape1_unique.txt", "w") as file:
    for lstfile in os.listdir(path_1):
        with open(path_1 + "/" + lstfile, "r") as openfile:
            for line in openfile:
                line = line.strip()
                line = line.split()
                unique1.extend(line)
    unique1 = list(set(unique1))
    for protein in unique1:
        file.writelines(protein + "\n")

# Read in result of GO Biological Process of size4_shape1 using reference list

size4shape1 = pd.read_csv("analysis_with_ref.csv")

size4shape1.columns = ["BP", "Ref_Num_Genes", "Motif_Num_Genes",
                       "Exp_Num_Genes", "Over_Or_Under", "Fold_Enrichment", "Raw_P_Val", "FDR_P_Val"]

size4shape1["Log_Fold_Enrichment"] = np.log(size4shape1["Fold_Enrichment"])

size4shape1["Log_FDR_P_Val"] = -np.log(size4shape1["FDR_P_Val"])
size4shape1["Log_P_Val_cmap"] = round(((size4shape1["Log_FDR_P_Val"] + np.min(
    size4shape1["Log_FDR_P_Val"])) * 255) / np.max(size4shape1["Log_FDR_P_Val"]))

plt.scatter(size4shape1["Log_Fold_Enrichment"],
            size4shape1["Log_FDR_P_Val"], s=size4shape1["Motif_Num_Genes"], ec="b", alpha=0.3)
# Add horizontal line to represent p-value = 1E-30
plt.axhline(30, c="r", linestyle="dashed")
plt.axvline(1, c="r", linestyle="dashed")

plt.title("GO Biological Process Enrichment in Size-4 Motif 1")
plt.xlabel("Log Fold Enrichment of Biological Process")
plt.ylabel("-Log FDR-Corrected p-Value")
plt.show()

plt.clf()

# Filter out highly significant and enriched proteins
size4shape1_filt = size4shape1[(size4shape1["Log_FDR_P_Val"] > 30) & (
    size4shape1["Log_Fold_Enrichment"] > 1)]
size4shape1_filt = size4shape1_filt.sort_values(by=["Motif_Num_Genes"])

plt.barh(y=size4shape1_filt["BP"], width=size4shape1_filt["Motif_Num_Genes"])
plt.ylabel("Biological Process")
plt.xlabel("Number of Proteins Enriched")
plt.title("Highly Enriched and Significant Proteins Found in Size-4 Motif 1")
plt.yticks(fontsize=6)
plt.show()
plt.clf()

# Concatenate files containing proteins in every motif occurrence for size-4 motif 2

path_2 = "./size4_shape2"

# Concatenate every row
unique2 = []
with open("size4_shape2_unique.txt", "w") as file:
    for lstfile in os.listdir(path_2):
        with open(path_2 + "/" + lstfile, "r") as openfile:
            for line in openfile:
                line = line.strip()
                line = line.split()
                unique2.extend(line)
    unique2 = list(set(unique2))
    for protein in unique2:
        file.writelines(protein + "\n")

# Read in result of GO Biological Process of size4_shape1 using reference list

size4shape2 = pd.read_csv("analysis_with_ref_2.csv")

size4shape2.columns = ["BP", "Ref_Num_Genes", "Motif_Num_Genes",
                       "Exp_Num_Genes", "Over_Or_Under", "Fold_Enrichment", "Raw_P_Val", "FDR_P_Val"]

size4shape2["Log_Fold_Enrichment"] = np.log(size4shape2["Fold_Enrichment"])

size4shape2["Log_FDR_P_Val"] = -np.log(size4shape2["FDR_P_Val"])
size4shape2["Log_P_Val_cmap"] = round(((size4shape2["Log_FDR_P_Val"] + np.min(
    size4shape2["Log_FDR_P_Val"])) * 255) / np.max(size4shape2["Log_FDR_P_Val"]))

plt.scatter(size4shape2["Log_Fold_Enrichment"],
            size4shape2["Log_FDR_P_Val"], s=size4shape2["Motif_Num_Genes"], ec="b", alpha=0.3)
# Add horizontal line to represent p-value = 1E-30
plt.axhline(20, c="r", linestyle="dashed")
plt.axvline(1, c="r", linestyle="dashed")

plt.title("GO Biological Process Enrichment in Size-4 Motif 2")
plt.xlabel("Log Fold Enrichment of Biological Process")
plt.ylabel("-Log FDR-Corrected p-Value")
plt.show()

plt.clf()

# Filter out highly significant and enriched proteins
size4shape2_filt = size4shape2[(size4shape2["Log_FDR_P_Val"] > 20) & (
    size4shape2["Log_Fold_Enrichment"] > 1)]
size4shape2_filt = size4shape2_filt.sort_values(by=["Motif_Num_Genes"])

plt.barh(y=size4shape2_filt["BP"], width=size4shape2_filt["Motif_Num_Genes"])
plt.ylabel("Biological Process")
plt.xlabel("Number of Proteins Enriched")
plt.title("Highly Enriched and Significant Proteins Found in Size-4 Motif 2")
plt.yticks(fontsize=6)
plt.show()
plt.clf()

# Read in result of PANTHER Pathway of size4_shape1 using reference list

size4shape1 = pd.read_csv("pathway_analysis_with_ref_1.csv")

size4shape1.columns = ["Pathway", "Ref_Num_Genes", "Motif_Num_Genes",
                       "Exp_Num_Genes", "Over_Or_Under", "Fold_Enrichment", "Raw_P_Val", "FDR_P_Val"]

size4shape1["Log_Fold_Enrichment"] = np.log(size4shape1["Fold_Enrichment"])

size4shape1["Log_FDR_P_Val"] = -np.log(size4shape1["FDR_P_Val"])


plt.scatter(size4shape1["Log_Fold_Enrichment"],
            size4shape1["Log_FDR_P_Val"], s=size4shape1["Motif_Num_Genes"], ec="b", alpha=0.3)
# Add horizontal line to represent p-value = 1E-30
plt.axhline(3, c="r", linestyle="dashed")
plt.axvline(1, c="r", linestyle="dashed")

plt.title("PANTHER Pathway Enrichment in Size-4 Motif 1")
plt.xlabel("Log Fold Enrichment of Pathway")
plt.ylabel("-Log FDR-Corrected p-Value")
plt.show()

plt.clf()

# Filter out highly significant and enriched proteins
size4shape1_filt = size4shape1[(size4shape1["Log_FDR_P_Val"] > 3) & (
    size4shape1["Log_Fold_Enrichment"] > 1)]

size4shape1_filt = size4shape1_filt.sort_values(by=["Motif_Num_Genes"])

plt.barh(y=size4shape1_filt["Pathway"],
         width=size4shape1_filt["Motif_Num_Genes"])
plt.ylabel("Pathway")
plt.xlabel("Number of Proteins Enriched")
plt.title("Highly Enriched and Significant Proteins Found in Size-4 Motif 1")
plt.yticks(fontsize=6)
plt.show()
plt.clf()


# Read in result of PANTHER Pathway of size4_shape2 using reference list

size4shape2 = pd.read_csv("pathway_analysis_with_ref_2.csv")

size4shape2.columns = ["Pathway", "Ref_Num_Genes", "Motif_Num_Genes",
                       "Exp_Num_Genes", "Over_Or_Under", "Fold_Enrichment", "Raw_P_Val", "FDR_P_Val"]

size4shape2["Log_Fold_Enrichment"] = np.log(size4shape2["Fold_Enrichment"])

size4shape2["Log_FDR_P_Val"] = -np.log(size4shape2["FDR_P_Val"])


plt.scatter(size4shape2["Log_Fold_Enrichment"],
            size4shape2["Log_FDR_P_Val"], s=size4shape2["Motif_Num_Genes"], ec="b", alpha=0.3)
# Add horizontal line to represent p-value = 1E-30
plt.axhline(3, c="r", linestyle="dashed")
plt.axvline(1, c="r", linestyle="dashed")

plt.title("PANTHER Pathway Enrichment in Size-4 Motif 2")
plt.xlabel("Log Fold Enrichment of Pathway")
plt.ylabel("-Log FDR-Corrected p-Value")
plt.show()

plt.clf()

# Filter out highly significant and enriched proteins
size4shape2_filt = size4shape2[(size4shape2["Log_FDR_P_Val"] > 3) & (
    size4shape2["Log_Fold_Enrichment"] > 1)]

size4shape2_filt = size4shape2_filt.sort_values(by=["Motif_Num_Genes"])

plt.barh(y=size4shape2_filt["Pathway"],
         width=size4shape2_filt["Motif_Num_Genes"])
plt.ylabel("Pathway")
plt.xlabel("Number of Proteins Enriched")
plt.title("Highly Enriched and Significant Proteins Found in Size-4 Motif 2")
plt.yticks(fontsize=6)
plt.show()
