from typing import List, Dict, Iterable
from collections import defaultdict
import sys

def de_bruijn_kmers(k_mers: List[str]) -> Dict[str, List[str]]:
    """Constructs a De Bruijn graph from a list of k-mers."""
    adj_list = defaultdict(list)
    for kmer in k_mers:
        node1, node2 = kmer[:-1], kmer[1:]
        adj_list[node1].append(node2)
    return dict(adj_list)

def calculate_degrees(graph: Dict[str, List[str]]) -> Dict[str, Dict[str, int]]:
    """Calculates the in-degree and out-degree of each node."""
    degrees = defaultdict(lambda: {'in': 0, 'out': 0})
    for node1, adj_nodes in graph.items():
        degrees[node1]['out'] += len(adj_nodes)
        for node2 in adj_nodes:
            degrees[node2]['in'] += 1
    return degrees

def contig_generation(Patterns: List[str]) -> List[str]:
    """Generates contigs from a set of k-mers."""
    graph = de_bruijn_kmers(Patterns)
    
    degrees = calculate_degrees(graph)
    
    contigs = []
    
    for node in graph:
        # If the node is a "branching" node (in-degree != 1 or out-degree != 1) or isolated
        if degrees[node]['in'] != 1 or degrees[node]['out'] != 1:
            if degrees[node]['out'] > 0:
                for neighbor in graph[node]:
                    contig = [node, neighbor]
                    # Follow the path until a branching or isolated node is reached
                    while degrees[neighbor]['in'] == 1 and degrees[neighbor]['out'] == 1:
                        next_node = graph[neighbor][0]
                        contig.append(next_node)
                        neighbor = next_node
                    contigs.append(''.join([contig[0]] + [c[-1] for c in contig[1:]]))
    
    return contigs