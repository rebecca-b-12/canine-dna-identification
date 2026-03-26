import argparse
import math
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from tabulate import tabulate

from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio.SeqRecord import SeqRecord

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for DNA identification program.

    Returns
    --------
    argparse.Namespace
        An object containing the parsed command-line arguments:
        - db (str): path to the FASTA database file
        - query (str): path to the query FASTA file
        - phylogeny (bool): Whether to display phylogenetic analysis
    """
    # Create the argument parser with a description of the program
    parser = argparse.ArgumentParser(
        description="Idnetify a DNA sequence by comparing it to a FASTA database")
    
    # Add argument for the FASTA database containing known sequences
    parser.add_argument(
        "--db", type=str, required=True, help="Path to FASTA database file (reference sequence)")
    

    # Add arguement for the query FASTA file to be identified
    parser.add_argument(
        "--query", type=str, required=True, help="Path to FATSA file containing the query sequence")

    # Optional flag: include phylogenetic analysis output
    parser.add_argument(
        "--phylogeny", action="store_true", help="Display phylogenetic tree (optional)")

    # Parse command-line arguments into a Namespace object
    args = parser.parse_args()

    # Check database and query file exist before proceeding
    if not os.path.exists(args.db):
        raise FileNotFoundError(f"Database file not found: {args.db}")
    if not os.path.exists(args.query):
        raise FileNotFoundError(f"Query file not found: {args.query}")
    
    # Return validated arguments
    return args

def sequence_identity(seq1: str, seq2: str) -> tuple[int, int]:
    """
    Calculate sequence identity between two DNA sequences.

    Parameters
    ----------
    seq1 : str
        First DNA sequence
    seq2 : str
        Secord DNA sequence

    Returns
    -------
    tuple[int, int]
        A tuple containing: 
        - matches (int): Number of matching nucleotide positions
        - valid_positions: Number of positions compared (excluding gaps or non-ATCG characters)
    
    """
    # Check sequences are not empty
    if not seq1 or not seq2:
        return 0, 0

    matches = 0 
    valid_positions = 0

    # Define valid DNA bases
    valid_bases = {"A", "T", "C", "G"}

    # Compare bases at each position in the sequences
    for base1, base2 in zip(seq1, seq2):

        # Only count positions where both bases are valid nucleotides
        if base1 in valid_bases and base2 in valid_bases:
            valid_positions += 1

            # Increase match counter if bases are identical
            if base1 == base2:
                matches += 1

    # Return total matches and number of valid comparisons
    return matches, valid_positions

def find_best_match(query_seq: str, database_records: list[SeqRecord]) -> tuple[SeqRecord, int, int, list[tuple[str, float, float]]]:
    """ 
    Identify the database sequence most similar to the query sequence.
    
    Parameters
    -----------
    query_seq : str
        DNA sequence to search for in the database.
    database_records : list[SeqRecord]
        List of reference sequences from the FASTA database.
    
    Returns
    -------
    tuple
        A tuple containing: 
        - best_record (SeqRecord): Most similar sequence in the database
        - best_match_count (int): Number of matching bases
        - best_valid_positions (int): Number of valid positions compared
        - similarity_results (list of tuples): 
            Each tuple contains (sequence_id, similarity, p_value)
    """

    # Check to ensure database is not empty
    if not database_records:
        raise ValueError("Database contains no sequences.")

    best_record = None
    best_similarity = -1
    best_match_count = 0
    best_valid_positions = 0

    similarity_results = []

    #Compare the query sequence with each sequence in the database
    for record in database_records:

        target_seq = str(record.seq)

        # Compute matches and valid positions
        matches, valid = sequence_identity(query_seq, target_seq)

        #Calculate similarity (proportion of matching bases)
        similarity = matches / valid if valid > 0 else 0

        # Calculate p-value
        p_value = math.exp(matches * math.log(0.25)) * len(database_records)

        # Store results for later use
        similarity_results.append((record.id, similarity, p_value))

        #Update best match if this sequence is more similar
        if similarity > best_similarity:
            best_similarity = similarity
            best_match_count = matches
            best_valid_positions= valid
            best_record = record

    return best_record, best_match_count, best_valid_positions, similarity_results

def extract_breed(description: str) -> str:
    """
    Extract the dog breed from a FASTA description string.

    Parameters
    ----------
    description : str
        Description line from a FASTA record.

    Returns
    -------
    str
        Extract breed name if present, otherwise "Unknown".
    """

    # Check for empty handle or missing description
    if not description:
        return "Unkown"

    # Search for breed annotation in square brackets
    match = re.search(r"\[breed=(.*?)\]", description)

    # Return extracted breed if found, otherwise default to "Unknown"
    return match.group(1) if match else "Unknown"

def build_tree(records: list[SeqRecord]) -> tuple[list[str], pd.DataFrame, DistanceMatrix]:
    """
    Construct a pairwise distance matrix and corresponding DistanceMatrix object.

    Parameters
    ----------
    records : list[SeqRecord]
        List of sequence records to compare
    
    Returns
    -------
    tuple
        - sequence_ids (list[str]): List of sequence identifiers
        - distance_df (pd.DataFrame): Full pairwise distance matrix
        - distance_matrix (DistanceMatrix): Lower-triangulat matrix for tree construction
    """

    # Check records are provided
    if not records: 
        raise ValueError("No records provided for tree construction.")

    sequence_ids = [r.id for r in records]

    # Initialise square DataFrame to store pairwise distances
    distance_df = pd.DataFrame(index=sequence_ids, columns=sequence_ids, dtype=float)

    # Compute pairwise distances (symmetrix matrix)
    for i, rec1 in enumerate(records):
        for j, rec2 in enumerate(records):

            # Avoid recomputing symmetric values
            if j < i:
                distance_df.iloc[i, j] = distance_df.iloc[j, i]
            
            seq1 = str(rec1.seq)
            seq2 = str(rec2.seq)

            matches, valid = sequence_identity(seq1, seq2)
            similarity = matches / valid if valid > 0 else 0

            distance = 1 - similarity

            distance_df.iloc[i, j] = distance
            distance_df.iloc[j, i] = distance

    # Convert full matrix into lower-triangular format 
    matrix = []
    for i in range(len(sequence_ids)):
        row = []
        for j in range(i + 1):
            row.append(float(distance_df.iloc[i, j]))
        matrix.append(row)

    # Create DistanceMatrix object for phylogenetic tree construction
    distance_matrix = DistanceMatrix(names=sequence_ids, matrix=matrix)

    return sequence_ids, distance_df, distance_matrix
    
def build_phylogeny(database_records: list[SeqRecord], best_id: str, probabilities: list[tuple[str, float, float]]) -> None:
    """
    Build and display phylogenetic trees from sequence data.

    Generates:
    - A full ASCII tree (all sequences) for terminal output
    - A filtered graphical tree (top 15 most similar sequences)

    Parameters
    ----------
    database_records : list[SeqRecord]
        All sequences from the database
    best_id : str
        ID of the best matching sequence
    probabilities : list[tuple[str, float, float]]
        List of (sequence_id, similarity, p-value) sorted by similarity
    """
    # Safety checks
    if not database_records:
        raise ValueError("No database records provided.")
    if best_id not in {r.id for r in database_records}:
        raise ValueError("Best match ID not found in database.")
    
    # Select top sequences
    top_n = 15

    # Create lookup to avoid repeated scans
    similarity_lookup = {seq_id: sim for seq_id, sim, _ in probabilities}

    top_ids = {seq_id for seq_id, _, _ in probabilities[:top_n]}
    top_ids.add(best_id)

    # Filter and sort records by similarity
    filtered_records = sorted([r for r in database_records if r.id in top_ids],
    key=lambda r: similarity_lookup[r.id],
    reverse=True)

    # Build full tree for terminal output
    all_ids, distance_df, distance_matrix = build_tree(database_records)

    print("\n" + "="*50)
    print("PHYLOGENETIC ANALYSIS")
    print("="*50)

    print(f"\nTotal sequences: {len(database_records)}")
    print("Method: Neighbour-Joining (distance = 1 - similarity)")
    print(f"\nBest match: {best_id}")

    print("\nDistance matrix (first 5x5 subset):\n")
    print(distance_df.iloc[:5, :5].round(4).to_string())

    constructor = DistanceTreeConstructor()
    full_tree = constructor.nj(distance_matrix)

    print("\n" + "-"*50)
    print("PHYLOGENETIC TREE")
    print("-"*50 + "\n")

    Phylo.draw_ascii(full_tree)

    # Build filtered tree for image output
    _, _, filtered_distance_matrix = build_tree(filtered_records)
    filtered_tree = constructor.nj(filtered_distance_matrix)

    # Highlight best match
    def label_func(clade):
        if clade.name == best_id:
            return f"{clade.name} ★"
        return clade.name

    # Plot tree
    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_subplot(1, 1, 1)

    Phylo.draw(
        filtered_tree,
        axes=ax,
        label_func=label_func,
        do_show=False,
        show_confidence=False
    )

    ax.set_xlabel("Genetic distance (1 - similarity)")
    ax.set_ylabel("Sequences (taxa)")

    # Save output
    output_file = "phylogenetic_tree.png"
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()

    print(f"\nTree image saved to: {output_file}")
    print("Note: Tree image shows top 15 sequences for clarity; full dataset used in analysis.")

def load_fasta_records(path: str) -> list[SeqRecord]:
    """
    Load DNA sequences from FASTA file.

    Parameters
    ----------
    path : str
        Path to the FASTA file
    
    Returns
    -------
    list[SeqRecord]
        List of parsed sequence records
    
    Raises
    ------
    FileNotFoundError
        If the file does not exist
    ValueError
        If the file contains no valid FASTA sequences
    """
    # Check the file exists before attempting to read it
    if not os.path.exists(path):
        raise FileNotFoundError("f:Fasta file '{path}' not found.")
    
    # Parse FASTA records using Biopython
    records = list(SeqIO.parse(path, "fasta"))

    # Ensure file is not empty or malformed
    if not records:
        raise ValueError(f"FASTA file '{path}' contains no sequences.")
    
    return records

def validate_query_sequence(query_record: SeqRecord) -> str:
    """
    Validate and extract a DNA sequence from a query record.

    Parameters
    ----------
    query_record : SeqRecord
        A sequence record containing the query DNA sequence
    
    Returns
    -------
    str
        Uppercase DNA sequence string
    
    Raises
    ------
    ValueError
        If the sequence is empty or contains no valid DNA bases (A, T, C, G)
    """
    
    # Convert sequence to uppercase string
    query_seq = str(query_record.seq.upper())

    # Check for empty sequence
    if not query_seq:
        raise ValueError("Query sequence is empty.")

    # Check that at least one valid DNA base is present
    if not any(base in "ATCG" for base in query_seq):
        raise ValueError("Query sequence contains no valid DNA bases (A, T, C, G).")
    
    return query_seq

def calculate_p_value(matches: int, database_size: int) -> float:
    """
    Calculate the probability (p-value) of observing a given number of matches
    by chance under a random DNA model.

    Assumes equal nucleotide probability (0.25 for A, T, C, G).

    Parameters
    ----------
    matches : int
        Number of matching nucleotide positions
    database_size : int
        Number of sequences in the database

    Returns
    -------
    float
        Estimated p-value for the observed similarity

    Raises
    ------
    ValueError
        If matches or database_size are negative
    """

    # Validate inputs
    if matches < 0:
        raise ValueError("Number of matches cannot be negative.")
    if database_size <= 0:
        raise ValueError("Database size must be positive.")

    # Compute probability under random model
    return (0.25 ** matches) * database_size

def print_best_match(best_record: SeqRecord, breed: str, difference: int) -> None:
    """
    Print details of the best matching sequence.

    Parameters
    ----------
    best_record : SeqRecord
        The sequence record with highest similarity to the query
    breed : str
        Extracted breed information from the sequence description
    difference : int
        Number of differing bases between query and best match
    """

    # Safety check
    if best_record is None:
        raise ValueError("Best match record is None.")

    print("\n" + "=" * 50)
    print("BEST MATCH")
    print("=" * 50)

    print(f"Closest sequence: {best_record.id}")
    print(f"Breed: {breed}")
    print(f"Number of differing bases: {difference}")

    # Convert sequence to string
    sequence = str(best_record.seq)

    # Handle short sequences safely
    preview = sequence[:50] + ("..." if len(sequence) > 50 else "")

    print(f"Sequence (first 50 bp): {preview}")

def print_similarity_table(probabilities: list[tuple[str, float, float]]) -> None:
    """
    Print a formatted table of sequence similarity results.

    Parameters
    ----------
    probabilities : list[tuple[str, float, float]]
        List of (sequence_id, similarity, p_value), sorted by similarity
    """

    # Safety check
    if not probabilities:
        print("\nNo similarity results to display.")
        return

    print("\n" + "=" * 50)
    print("SIMILARITY ACROSS ALL SEQUENCES")
    print("=" * 50)
    print("(Sorted by similarity to query sequence)\n")

    print(f"Total sequences: {len(probabilities)}\n")

    table_data = []

    for rank, (seq_id, similarity, p_value) in enumerate(probabilities, start=1):
        # Format p-value safely
        p_display = f"{p_value:.2e}" if p_value > 1e-300 else "<1e-300"

        # Highlight best match (rank 1)
        label = f"{seq_id} (best match)" if rank == 1 else seq_id

        table_data.append([
            rank,
            label,
            similarity,
            p_display])

    print(tabulate(
        table_data,
        headers=["Rank", "Sequence ID", "Similarity", "P-value"],
        tablefmt="fancy_grid",
        floatfmt=".4f"))

def print_summary(matches: int, valid: int, p_value: float) -> None:
    """
    Print a summary of the best match statistics.

    Parameters
    ----------
    matches : int
        Number of matching nucleotide positions
    valid : int
        Number of valid positions compared
    p_value : float
        Statistical significance of the match
    """

    # Safety checks
    if matches < 0 or valid < 0:
        raise ValueError("Matches and valid positions must be non-negative.")
    if valid == 0:
        raise ValueError("Valid positions cannot be zero when calculating similarity.")

    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)

    # Calculate similarity
    best_similarity = matches / valid

    # Format p-value safely
    p_display = f"{p_value:.2e}" if p_value > 1e-300 else "<1e-300"

    print(f"Best match similarity: {best_similarity:.4f}")
    print(f"Best match p-value: {p_display}")

def main() -> None:
    """
    Entry point for the DNA identification program.

    Workflow:
    1. Parse command-line arguments
    2. Load query and database sequences
    3. Validate query sequence
    4. Find best match in database
    5. Display results (best match, table, summary)
    6. Optionally generate phylogenetic tree
    """

    # Parse command-line arguments
    args = parse_arguments()

    # Load FASTA records
    query_records = load_fasta_records(args.query)
    database_records = load_fasta_records(args.db)

    # Ensure query file contains at least one sequence
    if not query_records:
        raise ValueError("Query FASTA file contains no sequences.")

    # Use the first query sequence
    query_record = query_records[0]

    # Validate query sequence
    query_sequence = validate_query_sequence(query_record)

    # Find best match
    best_record, matches, valid, similarity_results = find_best_match(
        query_sequence,
        database_records)

    if best_record is None:
        raise ValueError("No best match found.")

    # Extract additional info
    difference = valid - matches
    breed = extract_breed(best_record.description)
    p_value = calculate_p_value(matches, len(database_records))

    # Sort results by similarity (descending)
    similarity_results.sort(key=lambda x: x[1], reverse=True)

    # Output results
    print_best_match(best_record, breed, difference)
    print_similarity_table(similarity_results)
    print_summary(matches, valid, p_value)

    # Optional phylogenetic analysis
    if args.phylogeny:
        print(f"\nTotal sequences analysed: {len(database_records)}\n")

        build_phylogeny(
            database_records,
            best_record.id,
            similarity_results
        )


if __name__ == "__main__":
    main()