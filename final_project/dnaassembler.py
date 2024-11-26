from collections import defaultdict, deque
from typing import List, Dict, Set, Tuple, Optional
import random
import sys
import os
import math

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

def verify_eulerian_conditions(graph: Dict[str, List[str]]) -> Tuple[bool, Optional[str], Optional[str]]:
    """
    Verify if the graph has an Eulerian path or cycle.
    
    Args:
        graph: De Bruijn graph represented as adjacency list
        
    Returns:
        Tuple of (has_path, start_node, end_node)
    """
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    
    # Calculate degrees
    for node, neighbors in graph.items():
        out_degree[node] = len(neighbors)
        for neighbor in neighbors:
            in_degree[neighbor] += 1
    
    # Count nodes with unbalanced degrees
    start_nodes = []
    end_nodes = []
    for node in set(list(in_degree.keys()) + list(out_degree.keys())):
        diff = out_degree[node] - in_degree[node]
        if diff == 1:
            start_nodes.append(node)
        elif diff == -1:
            end_nodes.append(node)
        elif diff != 0:
            return False, None, None
    
    # Check if graph has Eulerian path
    if len(start_nodes) == 0 and len(end_nodes) == 0:
        # Eulerian cycle - any node can be start
        return True, next(iter(graph)), next(iter(graph))
    elif len(start_nodes) == 1 and len(end_nodes) == 1:
        # Eulerian path
        return True, start_nodes[0], end_nodes[0]
    
    return False, None, None

def build_de_bruijn_graph(kmers: Dict[str, int], check_eulerian: bool = True) -> Tuple[Dict[str, List[str]], bool]:
    """
    Build De Bruijn graph from k-mer counts.
    
    Args:
        kmers: Dictionary of k-mer counts
        check_eulerian: Whether to verify Eulerian conditions
        
    Returns:
        Tuple of (graph, is_eulerian)
    """
    edges = defaultdict(list)
    
    # Build graph with proper edge multiplicities
    for kmer, count in kmers.items():
        prefix, suffix = kmer[:-1], kmer[1:]
        edges[prefix].extend([suffix] * count)
    
    # Sort neighbors for deterministic results
    for node in edges:
        edges[node].sort()
    
    if check_eulerian:
        has_path, start, end = verify_eulerian_conditions(edges)
        if not has_path:
            return edges, False
    
    return edges, True

def find_eulerian_path(graph: Dict[str, List[str]]) -> List[str]:
    """
    Find Eulerian path in the graph using modified Hierholzer's algorithm.
    
    Args:
        graph: De Bruijn graph represented as adjacency list
        
    Returns:
        List of nodes forming Eulerian path
    """
    if not graph:
        return []
    
    # Verify Eulerian conditions and get proper start node
    has_path, start_node, _ = verify_eulerian_conditions(graph)
    if not has_path:
        return []
    
    # Copy graph to avoid modifying original
    graph = {k: v[:] for k, v in graph.items()}
    
    stack = [start_node]
    path = deque()
    
    while stack:
        current = stack[-1]
        if current in graph and graph[current]:
            next_node = graph[current].pop()
            stack.append(next_node)
        else:
            path.appendleft(stack.pop())
    
    return list(path)

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

def theoretical_reconstruction_possible(sequence_length: int, k: int, cyclic: bool = True) -> bool:
    """
    Determine if perfect reconstruction is theoretically possible.
    
    Args:
        sequence_length: Length of input sequence
        k: k-mer size
        cyclic: Whether sequence is circular
        
    Returns:
        Boolean indicating if reconstruction is possible
    """
    if not cyclic:
        return k <= sequence_length and (sequence_length - k + 1) >= sequence_length
    else:
        return k <= sequence_length and sequence_length >= k

def estimate_reconstruction_accuracy(sequence_length: int, k: int, alphabet_size: int = 4) -> float:
    """
    Estimate theoretical reconstruction accuracy.
    
    Args:
        sequence_length: Length of input sequence
        k: k-mer size
        alphabet_size: Size of alphabet (4 for DNA)
        
    Returns:
        Estimated accuracy between 0 and 1
    """
    if k <= 1:
        return 0.0
    
    possible_kmers = alphabet_size ** k
    actual_kmers = sequence_length
    
    if possible_kmers >= actual_kmers:
        coverage_ratio = actual_kmers / possible_kmers
        return min(1.0, (k - 1) / sequence_length * (1 - coverage_ratio))
    else:
        return possible_kmers / actual_kmers

def get_optimal_k(sequence_length: int, cyclic: bool = True) -> int:
    """
    Calculate optimal k-mer size for given sequence length.
    
    Args:
        sequence_length: Length of input sequence
        cyclic: Whether sequence is circular
        
    Returns:
        Optimal k-mer size
    """
    if cyclic:
        k = int(math.log(sequence_length, 4) + 1)
    else:
        k = int(math.log(sequence_length, 4) + 2)
    
    return max(2, min(k, sequence_length))

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

def reconstruct_from_kmers(sequence: str, k: int = None, cyclic: bool = True) -> Tuple[str, Dict]:
    """
    Main function to reconstruct sequence from k-mers.
    
    Args:
        sequence: Input sequence
        k: k-mer size (if None, uses optimal k)
        cyclic: Whether sequence is circular
        
    Returns:
        Tuple of (reconstructed sequence, details dictionary)
    """
    if k is None:
        k = get_optimal_k(len(sequence), cyclic)
    
    # Check if reconstruction is theoretically possible
    if not theoretical_reconstruction_possible(len(sequence), k, cyclic):
        estimated_accuracy = estimate_reconstruction_accuracy(len(sequence), k)
        return "", {
            "error": "Perfect reconstruction not theoretically possible",
            "estimated_accuracy": estimated_accuracy,
            "suggested_k": get_optimal_k(len(sequence), cyclic)
        }
    
    kmers = get_kmer_count_from_sequence(sequence, k, cyclic)
    graph, is_eulerian = build_de_bruijn_graph(kmers)
    
    if not is_eulerian:
        return "", {"error": "No Eulerian path exists in the graph"}
    
    path = find_eulerian_path(graph)
    reconstructed = reconstruct_sequence(path)
    
    score, details = score_sequence(sequence, reconstructed, circular=cyclic)
    details.update({
        'reconstructed_sequence': reconstructed,
        'estimated_accuracy': estimate_reconstruction_accuracy(len(sequence), k),
        'k_value': k,
        'sequence_length': len(sequence),
        'is_cyclic': cyclic
    })
    
    return reconstructed, details

def random_sequence(length: int, alphabet: str = 'ACGT') -> str:
    """Generate random sequence of given length"""
    return ''.join(random.choice(alphabet) for _ in range(length))

def main():
    """Main function for command line usage"""
    if len(sys.argv) < 2:
        print("Usage: python script.py <sequence_file> [k_value]")
        sys.exit(1)

    # Read sequence from file
    seq_file = sys.argv[1]
    with open(seq_file, "r") as handle:
        sequence = handle.read().strip().replace("\n", "")
        sequence = sequence.replace(" ", "")

    # Get k value from command line or use optimal
    k = int(sys.argv[2]) if len(sys.argv) > 2 else None
    
    print(f"Original sequence: {sequence[:100]}...")
    print(f"Length: {len(sequence)}")
    
    if k is None:
        k = get_optimal_k(len(sequence))
        print(f"Using optimal k = {k}")
    else:
        print(f"Using provided k = {k}")
    
    reconstructed, details = reconstruct_from_kmers(sequence, k=k, cyclic=True)
    
    if 'error' in details:
        print(f"\nError: {details['error']}")
        print(f"Estimated accuracy: {details['estimated_accuracy']:.1f}%")
        if 'suggested_k' in details:
            print(f"Suggested k value: {details['suggested_k']}")
    else:
        print(f"\nReconstructed sequence: {reconstructed[:100]}...")
        print(f"Score: {details['base_score']:.1f}")
        print(f"Identity: {details['percent_identity']:.1f}%")
        print(f"Rotation needed: {details.get('rotation', 0)} positions")
        print(f"Estimated accuracy: {details['estimated_accuracy']*100:.1f}%")

if __name__ == "__main__":
    main()