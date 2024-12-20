from collections import defaultdict, deque
from typing import List, Dict, Set, Tuple
import random
import sys
import os

def get_kmer_count_from_sequence(sequence: str, k: int = 3, cyclic: bool = True) -> Dict[str, int]:
    kmers = defaultdict(int)
    for i in range(len(sequence)):
        kmer = sequence[i:i + k]
        if len(kmer) != k:
            if cyclic:
                kmer += sequence[:(k - len(kmer))]
            else:
                continue
        kmers[kmer] += 1
    return dict(kmers)

def build_de_bruijn_graph(kmers: Dict[str, int]) -> Dict[str, List[str]]:
    edges = defaultdict(list)
    for kmer, count in kmers.items():
        prefix, suffix = kmer[:-1], kmer[1:]
        edges[prefix].extend([suffix] * count)
    return edges

def find_eulerian_path(graph: Dict[str, List[str]]) -> List[str]:
    if not graph:
        return []
        
    # Calculate in/out degrees
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    for node, neighbors in graph.items():
        out_degree[node] = len(neighbors)
        for neighbor in neighbors:
            in_degree[neighbor] += 1
            
    # Find start node
    start_node = next((node for node in out_degree if out_degree[node] > in_degree[node]), 
                     next(iter(graph)))
    
    # Find path using Hierholzer's algorithm
    stack = [start_node]
    path = deque()
    
    while stack:
        node = stack[-1]
        if graph[node]:
            next_node = graph[node].pop()
            stack.append(next_node)
        else:
            path.appendleft(stack.pop())
            
    return list(path)

def reconstruct_sequence(path: List[str]) -> str:
    if not path:
        return ""
    return path[0] + ''.join(node[-1] for node in path[1:])

def score_sequence(original: str, reconstructed: str, circular: bool = True,
                  match_score: float = 1.0, mismatch_score: float = -0.5,
                  length_penalty: float = -1.0) -> Tuple[float, Dict]:
    if not circular:
        min_length = min(len(original), len(reconstructed))
        length_diff = abs(len(original) - len(reconstructed))
        matches = sum(1 for i in range(min_length) if original[i] == reconstructed[i])
        mismatches = min_length - matches
        
        base_score = (matches * match_score) + (mismatches * mismatch_score)
        length_penalty_score = length_diff * length_penalty
        
        details = {
            'matches': matches,
            'mismatches': mismatches,
            'length_diff': length_diff,
            'base_score': base_score,
            'length_penalty': length_penalty_score,
            'percent_identity': (matches / min_length * 100) if min_length > 0 else 0,
            'mode': 'linear'
        }
        return base_score + length_penalty_score, details
    
    # For circular sequences, try all rotations
    best_score = float('-inf')
    best_details = {}
    
    for i in range(len(reconstructed)):
        rotated = reconstructed[i:] + reconstructed[:i]
        score, details = score_sequence(original, rotated, False, 
                                     match_score, mismatch_score, length_penalty)
        if score > best_score:
            best_score = score
            best_details = details
            best_details.update({
                'rotation': i,
                'mode': 'circular',
                'aligned_sequence': rotated
            })
    
    return best_score, best_details

def reconstruct_from_kmers(sequence: str, k: int, cyclic: bool = True) -> Tuple[str, Dict]:
    kmers = get_kmer_count_from_sequence(sequence, k, cyclic)
    graph = build_de_bruijn_graph(kmers)
    path = find_eulerian_path(graph)
    reconstructed = reconstruct_sequence(path)
    
    score, details = score_sequence(sequence, reconstructed, circular=cyclic)
    details['reconstructed_sequence'] = reconstructed
    return reconstructed, details

def random_sequence(length: int) -> str:
    return ''.join(random.choice('ACGT') for _ in range(length))

if __name__ == "__main__":
    # Example usage
    # read seq from command line (it will be a txt file)
    # read the file and save the text as a string
    seq_file = sys.argv[1]
    with open(seq_file, "r") as handle:
        sequence = handle.read().strip()

    # make sure the sequence is a string with one line (no new line characters)
    sequence = sequence.replace("\n", "")
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\t", "")
    sequence = sequence.replace("\r", "")


    # sequence = random_sequence(10000)
    print(f"Original sequence: {sequence[:100]}...")
    print(f"Length: {len(sequence)}")

    kval = int(float(sys.argv[2]))
    
    reconstructed, details = reconstruct_from_kmers(sequence, k=kval, cyclic=True)
    print(f"K value: {kval}")
    print(f"\nReconstructed sequence: {reconstructed[:100]}...")
    print(f"Score: {details['base_score']:.1f}")
    print(f"Identity: {details['percent_identity']:.1f}%")
    print(f"Rotation needed: {details.get('rotation', 0)} positions")