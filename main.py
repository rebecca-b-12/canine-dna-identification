import argparse
import re
import pandas as pd
import math

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
from tabulate import tabulate

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for DNA identification program.

    Returns
    argparse.Namespace
        An object containing the parsed command-line arguments:
        - db: path to the FASTA database file
        - query: path to the query FASTA file
    """
    #Create the argument parser with a description of the program
    parser = argparse.ArgumentParser(description="DNA Identification program")
    #Add argument for the FASTA database containing known sequences
    parser.add_argument("--db", required=True, help="FASTA database file")
    #Add arguement for the query FASTA file to be identified
    parser.add_argument("--query", required=True, help="Query FASTA file")

    parser.add_argument("--phylogeny", action="store_true", help="Display phylogenetic tree")
    #Parse and return the command-line arguments
    return parser.parse_args()

def sequence_identity(seq1: str, seq2: str) -> tuple[int, int]:
    """
    Compare two DNA sequences and count matching bases.

    Parameters
    seq1 : str
        First DNA sequence
    seq2 : str
        Secord DNA sequence

    Returns
    tuple[int, int]
        A tuple containing: 
        - matches: number of matching nucleotide positions
        - valid_positions: number of positions compared (excluding gaps or ambiguous bases)
    
    """

    matches = 0 
    valid_positions = 0

    #Compare bases at each position in the sequences
    for a, b in zip(seq1, seq2):
        #Only count positions where both bases are valid nucleotides
        if a in "ATCG" and b in "ATCG":
            valid_positions += 1
            #Increase match counter if bases are identical
            if a == b:
                matches += 1

    return matches, valid_positions

def find_best_match(query_seq: str, database_records: list[SeqRecord]) -> tuple[SeqRecord, int, int, list[tuple[str, float, float]]]:
    """ 
    Identify the database sequence most similar to the query sequence.
    
    Parameters
    query_seq : str
        DNA sequence to search for in the database.
    database_records : list[SeqRecord]
        List of sequence records from the FASTA database.
    
    Returns
    tuple
        A tuple containing: 
        - best_record: the SeqRecord with the highest similarity
        - best_score: number of matching bases
        - best_valid: number of valid positions compared
        - probabilities: list of (sequence_id, similarity_probability)
    """
    best_record = None
    best_similarity = -1
    best_matches = 0
    best_valid = 0

    probabilities = []

    #Compare the query sequence with each sequence in the database
    for record in database_records:

        sequence = str(record.seq)
        matches, valid = sequence_identity(query_seq, sequence)

        #Calculate similarity probability
        similarity = matches / valid if valid > 0 else 0
        p_value = math.exp(matches * math.log(0.25)) * len(database_records)
        probabilities.append((record.id, similarity, p_value))

        #Update best match if this sequence has a higher score
        if similarity > best_similarity:
            best_similarity = similarity
            best_matches = matches
            best_valid = valid
            best_record = record

    return best_record, best_matches, best_valid, probabilities

def extract_breed(description: str) -> str:
    """
    Extract the dog breed from a FASTA description line.

    Parameters
    description : str
        Description string from a FASTA record.

    Returns
    str
        The breed name if found, otherwise "Unknown".
    """

    #Search for pattern like "[breed=boxer]"
    match = re.search(r"\[breed=(.*?)\]", description)
    if match:
        return match.group(1)
    else:
        return "Unknown"
    
def build_phylogeny(database_records: list[SeqRecord], best_id: str) -> None:

    ids = [record.id for record in database_records]

    df = pd.DataFrame(index=ids, columns=ids)

    for i, rec1 in enumerate(database_records):
        for j, rec2 in enumerate(database_records):

            seq1 = str(rec1.seq)
            seq2 = str(rec2.seq)

            matches, valid = sequence_identity(seq1, seq2)

            similarity = matches / valid if valid > 0 else 0
            distance = 1 - similarity

            df.iloc[i, j] = distance
    
    print("\nDistance matrix (used to construct phylogenetic tree):")
    print("(Distance = 1 - similarity; showing first 5 rows for readability)\n")
    print(df.head())
    print(f"\nMatrix size: {df.shape[0]} x {df.shape[1]}")
    print("These distances are used to build the Neighbour-Joining tree.\n")

    print("\nPhylogenetic tree (Neighbour-Joining):")
    print("All sequences are included.")
    print("Branch lengths represent sequence dissimilarity (1 - similarity).")
    print(f"Best match: {best_id}\n")

    matrix = []
    for i in range(len(ids)):
        row = []
        for j in range(i + 1):
            row.append(float(df.iloc[i, j]))
        matrix.append(row)
    
    distance_matrix = DistanceMatrix(names=ids, matrix=matrix)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    Phylo.draw_ascii(tree)

def load_fasta_records(path: str) -> list[SeqRecord]:
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"FASTA file '{path}' contains no sequences.")
    return records

def validate_query_sequence(query_record: SeqRecord) -> str:
    query_seq = str(query_record.seq.upper())

    if not any(base in "ATCG" for base in query_seq):
        raise ValueError("Query sequence contains no valid DNA bases (A, T, C, G).")
    
    return query_seq

def calculate_p_value(matches: int, db_size: int) -> float:
    return (0.25 ** matches) * db_size

def print_best_match(best_record: SeqRecord, breed: str, difference: int) -> None:
    print("\n" + "="*50)
    print("BEST MATCH")
    print("="*50)

    print("Closest sequence:", best_record.id)
    print("Breed:", breed)
    print("Difference (number of differing bases):", difference)

    seq = str(best_record.seq)
    print("Sequence (first 50 bp):", seq[:50] + "...")

def print_similarity_table(probabilities: list[tuple[str, float]]) -> None:
    print("\n" + "="*50)
    print("SIMILARITY ACROSS ALL SEQUENCES")
    print("="*50)
    print("(Sorted by similarity to query sequence)\n")

    print(f"Total sequences: {len(probabilities)}\n")

    table_data = []

    for i, (seq_id, similarity, p_value) in enumerate(probabilities, start=1):
        p_display = f"{p_value:.2e}" if p_value > 1e-300 else "<1e-300"
        table_data.append([i, seq_id, similarity, p_display])

    print(tabulate(
        table_data, 
        headers=["Rank", "Sequence ID", "Similarity", "P-value"], 
        tablefmt="fancy_grid",
        floatfmt=".4f"))

def print_summary(matches: int, valid: int, p_value: float) -> None:
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)

    best_similarity = matches / valid if valid > 0 else 0
    p_display = f"{p_value:.2e}" if p_value > 1e-300 else "<1e-300"

    print(f"Best match similarity: {best_similarity:.4f}")
    print(f"Best match p-value: {p_display}")

def main() -> None:
    args = parse_arguments()

    query_records = load_fasta_records(args.query)
    database_records = load_fasta_records(args.db)

    query_record = query_records[0]

    query_seq = validate_query_sequence(query_record)

    best_record, matches, valid, probabilities = find_best_match(query_seq, database_records)

    difference = valid - matches
    breed = extract_breed(best_record.description)
    p_value = calculate_p_value(matches, len(database_records))

    probabilities.sort(key=lambda x: x[1], reverse=True)

    print_best_match(best_record, breed, difference)
    print_similarity_table(probabilities)
    print_summary(matches, valid, p_value)

    if args.phylogeny:
        print("\n" + "="*50)
        print("PHYLOGENETIC TREE")
        print("="*50)

        print(f"\nTotal sequences analysed: {len(database_records)}\n")

        build_phylogeny(database_records, best_record.id)

if __name__ == "__main__":
    main() 