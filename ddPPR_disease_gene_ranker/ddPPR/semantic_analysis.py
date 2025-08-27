import networkx as nx
import numpy as np
import math
import pandas as pd
from itertools import combinations
from collections import defaultdict
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

def compute_information_content(G):
    """Compute IC for each node based on descendant count (inverse specificity)."""
    N = len(G)
    IC = {}
    for node in G.nodes:
        freq = len(nx.descendants(G, node)) + 1  # +1 to avoid log(0); retrieves how many nodes are accessible from a given node
        IC[node] = -math.log(freq / N) # info content is the log of how many nodes are reachable from that node
    return IC

def get_LCAs(G, node1, node2):
    """Return the set of lowest common ancestors of two nodes in a DAG."""
    anc1 = nx.ancestors(G, node1) | {node1}
    anc2 = nx.ancestors(G, node2) | {node2}
    return anc1 & anc2

def resnik_similarity(G, IC, node1, node2):
    """Compute Resnik similarity between two nodes."""
    LCAs = get_LCAs(G, node1, node2)
    if not LCAs:
        return 0.0
    # Get the highest LCA (the node that can reach the farthest in the network, the highest node in the hirearchy)
    return max(IC[a] for a in LCAs) # IC being a dictionary for all nodes in the graph

def compute_resnik_matrix(G, pathways, IC, normalize=True):
    """Compute a symmetric Resnik similarity matrix between pathways."""
    max_IC = max(IC.values()) if normalize else 1.0
    matrix = pd.DataFrame(index=pathways, columns=pathways, dtype=float)

    for a in pathways:
        for b in pathways:
            if pd.isna(matrix.loc[a, b]):
                sim = resnik_similarity(G, IC, a, b)
                if normalize and max_IC > 0:
                    sim /= max_IC
                matrix.loc[a, b] = sim
                matrix.loc[b, a] = sim  # symmetry
# Fill diagonal with self-similarity (should be highest possible value)
    for a in pathways:
        matrix.loc[a, a] = 1.0  # because sim(node, node) = max_IC / max_IC = 1.0 if normalized

    return matrix

def filter_redundant_pathways_by_resnik(G, enriched_pathways, similarity_threshold=0.3):
    """
    Given a DAG G and a list of enriched pathways (nodes in G),
    return a filtered list keeping only the most specific pathway per cluster
    based on Resnik similarity.
    
    Args:
        G (networkx.DiGraph): Reactome DAG.
        enriched_pathways (list): List of Reactome pathway IDs (R-HSA-XXXXX).
        similarity_threshold (float): Distance threshold for clustering (1 - similarity).
    
    Returns:
        filtered_pathways (list): Most specific (high IC) pathway per cluster.
        resnik_matrix (pd.DataFrame): Full similarity matrix (normalized).
    """
    # Step 1: Compute information content
    IC = compute_information_content(G)

    # Step 2: Compute Resnik similarity matrix
    resnik_matrix = compute_resnik_matrix(G=G, pathways=enriched_pathways, IC=IC, normalize=True)

    # Step 3: Cluster similar pathways using hierarchical clustering
    dist_matrix = 1 - resnik_matrix.fillna(0)
    np.fill_diagonal(dist_matrix.values, 0)  # Fix the diagonal
    condensed = squareform(dist_matrix.values)
    linkage_matrix = linkage(condensed, method='average') # the clustering process itself

    cluster_ids = fcluster(linkage_matrix, t=similarity_threshold, criterion='distance')
    # Step 4: For each cluster, keep the most specific (highest IC) term
    clusters = defaultdict(list)
    for pathway, cluster_id in zip(enriched_pathways, cluster_ids):
        clusters[cluster_id].append(pathway)

    filtered_pathways = []
    for group in clusters.values():
        most_specific = max(group, key=lambda p: IC[p])
        filtered_pathways.append(most_specific)

    return filtered_pathways, resnik_matrix