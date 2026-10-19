import pandas as pd
from utils.selection_utils import ordered_edges


def select_top_edges_per_seed(edges, k_edges, edge_random_generator):
    edges = ordered_edges(edges)
    if k_edges is None or len(edges) <= k_edges:
        return edges
    boundary = edges.iloc[k_edges - 1]["weight"]
    edges_above_boundary = edges[edges["weight"] > boundary]
    edges_on_boundary = edges[edges["weight"] == boundary] 
    slots_for_edges_on_boundary = k_edges - len(edges_above_boundary)
    if slots_for_edges_on_boundary == len(edges_on_boundary):
        return edges.head(k_edges).reset_index(drop=True)
    random_order = edge_random_generator.permutation(len(edges_on_boundary))
    chosen_positions = random_order[:slots_for_edges_on_boundary]
    chosen_boundary_edges = edges_on_boundary.iloc[chosen_positions] # draw remaining slots from edges on the boundray with equal weights
    selected_edges = pd.concat([edges_above_boundary, chosen_boundary_edges], ignore_index=True)
    return ordered_edges(selected_edges)
