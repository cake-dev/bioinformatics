from collections import defaultdict, deque
from typing import List, Dict, Set, Tuple, Optional
import random
import sys
import os
import math
from statistics import mean, median

def get_kmer_count_from_sequence(sequence: str, k: int = 3, cyclic: bool = True) -> Dict[str, int]:
    """
    Generate k-mer counts from a sequence.
    
    Args:
        sequence: Input DNA sequence
        k: k-mer size
        cyclic: Whether to treat sequence as circular
        
    Returns:
        Dictionary of k-mer counts
    """
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

def reconstruct_sequence(path: List[str]) -> str:
    """
    Reconstruct sequence from path in De Bruijn graph.
    
    Args:
        path: List of nodes forming path
        
    Returns:
        Reconstructed sequence
    """
    if not path:
        return ""
    return path[0] + ''.join(node[-1] for node in path[1:])

def score_sequence(original: str, reconstructed: str, circular: bool = True,
                  match_score: float = 1.0, mismatch_score: float = -0.5,
                  length_penalty: float = -1.0) -> Tuple[float, Dict]:
    """
    Score reconstruction against original sequence.
    
    Args:
        original: Original sequence
        reconstructed: Reconstructed sequence
        circular: Whether sequence is circular
        match_score: Score for matching bases
        mismatch_score: Penalty for mismatching bases
        length_penalty: Penalty for length differences
        
    Returns:
        Tuple of (score, details dictionary)
    """
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

def analyze_sequence_complexity(sequence: str, window_size: int = 100) -> float:
    """
    Analyze sequence complexity to help determine optimal parameters.
    Returns a complexity score between 0 and 1.
    """
    complexity_scores = []
    for i in range(0, len(sequence) - window_size, window_size):
        window = sequence[i:i + window_size]
        unique_kmers = len(set(window[j:j+3] for j in range(len(window)-2)))
        max_unique = min(window_size-2, 64)  # 64 is max possible 3-mers
        complexity_scores.append(unique_kmers / max_unique)
    return mean(complexity_scores) if complexity_scores else 0.5

def get_optimal_params(sequence: str) -> Dict:
    """
    Determine optimal parameters based on sequence characteristics.
    """
    length = len(sequence)
    complexity = analyze_sequence_complexity(sequence)
    
    # Base k calculation using sequence length
    base_k = int(math.log(length, 4)) + 1
    
    # Adjust k based on sequence complexity
    if complexity < 0.3:  # Very repetitive
        k = max(5, base_k - 2)
    elif complexity < 0.7:  # Moderate complexity
        k = base_k
    else:  # High complexity
        k = base_k + 1
        
    # Ensure k is within reasonable bounds
    k = max(4, min(k, int(math.log2(length))))
    
    # Calculate filtering threshold based on complexity
    noise_threshold = 0.05 + (0.15 * (1 - complexity))
    
    return {
        'k': k,
        'noise_threshold': noise_threshold,
        'sequence_complexity': complexity
    }

def build_weighted_de_bruijn_graph(kmers: Dict[str, int]) -> Tuple[Dict[str, List[Tuple[str, int]]], bool]:
    """
    Build De Bruijn graph with edge weights based on k-mer frequencies.
    """
    edges = defaultdict(list)
    frequencies = defaultdict(int)
    
    # Calculate average and median frequencies
    avg_freq = mean(kmers.values())
    med_freq = median(kmers.values())
    threshold = min(avg_freq, med_freq) * 0.1  # Adaptive threshold
    
    # Build graph with frequencies
    for kmer, count in kmers.items():
        if count >= threshold:
            prefix, suffix = kmer[:-1], kmer[1:]
            edges[prefix].append((suffix, count))
            frequencies[prefix] += count
            frequencies[suffix] += count
    
    # Sort edges by weight for each node
    for node in edges:
        edges[node].sort(key=lambda x: x[1], reverse=True)
    
    # Check if graph is approximately Eulerian
    has_path = check_approximate_eulerian(edges, frequencies)
    
    return edges, has_path

def check_approximate_eulerian(graph: Dict[str, List[Tuple[str, int]]], 
                             frequencies: Dict[str, int]) -> bool:
    """
    Check if graph is approximately Eulerian (allows small imbalances).
    """
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    for node, edges in graph.items():
        out_sum = sum(weight for _, weight in edges)
        out_degree[node] = out_sum
        for next_node, weight in edges:
            in_degree[next_node] += weight
    
    # Allow small imbalances based on frequency
    imbalanced_nodes = 0
    for node in set(in_degree.keys()) | set(out_degree.keys()):
        diff = abs(out_degree[node] - in_degree[node])
        if diff > 0.1 * frequencies[node]:  # 10% tolerance
            imbalanced_nodes += 1
    
    return imbalanced_nodes <= 2  # Allow up to 2 imbalanced nodes

def find_best_path(graph: Dict[str, List[Tuple[str, int]]], 
                   start_node: Optional[str] = None) -> List[str]:
    """
    Find highest-weight path through the graph.
    """
    if not graph:
        return []
    
    # Copy graph
    graph = {k: v[:] for k, v in graph.items()}
    
    # Find start node with highest outgoing weight if not provided
    if not start_node:
        start_node = max(graph.keys(), 
                        key=lambda x: sum(w for _, w in graph.get(x, [])))
    
    path = []
    current = start_node
    
    while current in graph and graph[current]:
        path.append(current)
        # Take highest weight edge
        next_node, _ = max(graph[current], key=lambda x: x[1])
        # Remove used edge
        graph[current] = [(n, w) for n, w in graph[current] if n != next_node]
        current = next_node
    
    path.append(current)
    return path

def reconstruct_from_kmers(sequence: str, k: int = None) -> Tuple[str, Dict]:
    """
    Improved sequence reconstruction with automatic parameter optimization.
    """
    # Get optimal parameters if not specified
    if k is None:
        params = get_optimal_params(sequence)
        k = params['k']
    else:
        params = {'k': k, 'sequence_complexity': analyze_sequence_complexity(sequence)}
    
    # Get k-mer counts
    kmers = get_kmer_count_from_sequence(sequence, k, cyclic=True)
    
    # Build weighted graph
    graph, is_eulerian = build_weighted_de_bruijn_graph(kmers)
    
    # Find best path
    path = find_best_path(graph)
    reconstructed = reconstruct_sequence(path)
    
    # Score reconstruction
    score, details = score_sequence(sequence, reconstructed, circular=True)
    
    # Add analysis details
    details.update({
        'reconstructed_sequence': reconstructed,
        'k_value': k,
        'sequence_length': len(sequence),
        'sequence_complexity': params['sequence_complexity'],
        'is_cyclic': True
    })
    
    return reconstructed, details

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <sequence_file> [k_value]")
        sys.exit(1)

    # Read sequence
    seq_file = sys.argv[1]
    with open(seq_file, "r") as handle:
        sequence = handle.read().strip().replace("\n", "")
        sequence = sequence.replace(" ", "")
        sequence = sequence.replace("\t", "")

    with open("original_seq.txt", "w") as handle2:
        handle2.write(sequence)

    # Get k value
    k = int(sys.argv[2]) if len(sys.argv) > 2 else None
    
    print(f"Original sequence: {sequence[:100]}...")
    print(f"Length: {len(sequence)}")
    
    # Get sequence complexity
    complexity = analyze_sequence_complexity(sequence)
    print(f"Sequence complexity score: {complexity:.2f}")
    
    if k is None:
        params = get_optimal_params(sequence)
        k = params['k']
        print(f"Using optimal k = {k} (based on sequence characteristics)")
    else:
        print(f"Using provided k = {k}")
    
    reconstructed, details = reconstruct_from_kmers(sequence, k=k)
    
    print(f"\nReconstructed sequence: {reconstructed[:100]}...")
    print(f"Score: {details['base_score']:.1f}")
    print(f"Identity: {details['percent_identity']:.1f}%")
    print(f"Rotation needed: {details.get('rotation', 0)} positions")

    # save both the original and the reconstructed sequence to files
    with open("original_seq.txt", "w") as handle:
        handle.write(sequence)
    with open("reconstructed_seq.txt", "w") as handle:
        handle.write(reconstructed)

if __name__ == "__main__":
    main()