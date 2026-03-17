import argparse
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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

def find_best_match(query_seq: str, database_records: list[SeqRecord]) -> tuple[SeqRecord, int, int, list[tuple[str, float]]]:
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
    best_score = -1
    best_valid = 0
    probabilities = []

    #Compare the query sequence with each sequence in the database
    for record in database_records:

        sequence = str(record.seq)
        matches, valid = sequence_identity(query_seq, sequence)

        #Calculate similarity probability
        probability = matches / valid if valid > 0 else 0
        probabilities.append((record.id, probability))

        #Update best match if this sequence has a higher score
        if matches > best_score:
            best_score = matches
            best_valid = valid
            best_record = record

    return best_record, best_score, best_valid, probabilities

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


def main() -> None:
    """
    Run the DNA identification workflow

    The program:
    1. Parses command-line arguments.
    2. Loads the query sequence and database sequences.
    3. Identifies the closest matching sequence.
    4. Calculates sequence differences and probabilities.
    5. Prints the results.
    """

    #Parse command-line inputs
    args = parse_arguments()

    #Load FASTA records
    query_records = list(SeqIO.parse(args.query, "fasta"))

    #Check if query record is empty, raise error if empty.
    if not query_records:
        raise ValueError("Query FASTA file contains no sequences.")
    query_record = query_records[0]

    database_records = list(SeqIO.parse(args.db, "fasta"))

    #Check if database is empty, raise error if empty.
    if not database_records:
        raise ValueError("Database FASTA file contains no sequences.")

    #convert query sequence into string for comparison
    query_seq = str(query_record.seq).upper()
    
    #Ensure query sequence contains valid DNA bases
    if not any(base in "ATCG" for base in query_seq):
        raise ValueError("Query sequence contains no valid DNA bases (A, T, C, G).")

    #Identify the closest sequence match
    best_record, matches, valid, probabilities = find_best_match(query_seq, database_records)
    
    difference = valid - matches

    #Extract breed information from FASTA description
    breed = extract_breed(best_record.description)

    #Estimate p-value for observing this similarity by chance
    p_value = (0.25 ** matches) * len(database_records)

    print("Closest sequence:", best_record.id)
    print("Breed:", breed)
    print("Difference:", difference)

    probabilities.sort(key=lambda x: x[1], reverse=True)

    print("\nProbabilities across database:")

    max_id_length = max(len(seq_id) for seq_id, _ in probabilities)
    print(f"{'Rank':<5}{'Sequence ID':<{max_id_length}} Probability")
    print("-" * (max_id_length + 20))

    for i, (seq_id, probability) in enumerate(probabilities, start=1):
        print(f"{i:>2}. {seq_id:<20} probability={probability:.4f}")
    
    print("\np_value:", p_value)



if __name__ == "__main__":
    main() 