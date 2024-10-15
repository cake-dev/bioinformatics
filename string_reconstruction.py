from typing import List, Dict, Iterable
import sys
from collections import defaultdict

# Constructs a genome from a path of k-mers
def genome_path(path: List[str]) -> str:
    """Reconstructs a genome from a given path of k-mers."""
    genome = path[0]
    for pattern in path[1:]:
        genome += pattern[-1]
    return genome

# Constructs a De Bruijn graph from a list of k-mers
def de_bruijn_kmers(k_mers: List[str]) -> Dict[str, List[str]]:
    """Constructs a De Bruijn graph from a list of k-mers."""
    adj_list = defaultdict(list)
    for kmer in k_mers:
        node1, node2 = kmer[:-1], kmer[1:]
        adj_list[node1].append(node2)
    return dict(adj_list)

# Extends the Eulerian cycle in the graph
def extend_cycle(cycle: List[str], marked_graph: Dict[str, List[str]]) -> List[str]:
    """Extends the Eulerian cycle from a given node in the marked graph."""
    if cycle:
        cycle.pop()  # remove the repeated node at the end
        new_start_index = next(i for i, node in enumerate(cycle) if node in marked_graph)
        cycle = cycle[new_start_index:] + cycle[:new_start_index]
        cycle.append(cycle[0])  # re-add the repeated node
        current_node = cycle[-1]
    else:
        current_node = next(iter(marked_graph))  # get an arbitrary node from the graph
        cycle = [current_node]
    
    while current_node in marked_graph:
        old_node = current_node
        current_node = marked_graph[old_node].pop()
        if not marked_graph[old_node]:
            del marked_graph[old_node]  # remove the node if no more edges
        cycle.append(current_node)
    
    return cycle

# Constructs an Eulerian cycle in a graph
def eulerian_cycle(g: Dict[str, List[str]]) -> List[str]:
    """Constructs an Eulerian cycle in a graph. Assumes the graph is Eulerian and connected."""
    cycle = []
    while g:
        cycle = extend_cycle(cycle, g)
    return cycle

# Fixes unbalanced nodes in a graph to make it Eulerian
def fix_unbalanced(g: Dict[str, List[str]]) -> tuple[str, str]:
    """Finds and fixes unbalanced nodes in the graph."""
    total_degree = defaultdict(int)
    
    for node1, adj_nodes in g.items():
        for node2 in adj_nodes:
            total_degree[node1] += 1  # Out-degree
            total_degree[node2] -= 1  # In-degree

    s, t = None, None
    for node, tot_degree in total_degree.items():
        if tot_degree == 1:
            t = node
        elif tot_degree== -1:
            s = node

    if s and t:
        g.setdefault(s, []).append(t)
    
    return s, t

# Constructs an Eulerian path in a graph
def eulerian_path(g: Dict[str, List[str]]) -> List[str]:
    """Constructs an Eulerian path in a graph, assuming the graph is nearly Eulerian."""
    s, t = fix_unbalanced(g)
    cycle = eulerian_cycle(g)
    
    if s:
        cycle.pop()  # Remove the duplicate last node
        t_index = next(i for i, (u, v) in enumerate(zip(cycle, cycle[1:])) if u == s and v == t)
        cycle = cycle[t_index + 1:] + cycle[:t_index + 1]
    
    return cycle

# Reconstructs a string from its k-mer composition
def string_reconstruction(patterns: List[str], k: int) -> str:
    """Reconstructs a string from its k-mer composition."""
    pass
