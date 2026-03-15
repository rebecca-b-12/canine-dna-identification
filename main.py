import argparse
import re
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="DNA Identification Service")
    parser.add_argument("--db", required=True, help="FASTA database file")
    parser.add_argument("--query", required=True, help="Query FASTA file")
    return parser.parse_args()

def sequence_identity(seq1, seq2):
    return sum(a == b for a, b in zip(seq1, seq2))

def find_best_match(query_seq, database_records):

    best_record = None
    best_score = -1

    for record in database_records:

        sequence = str(record.seq)
        score = sequence_identity(query_seq, sequence)

        if score > best_score:
            best_score = score
            best_record = record
    return best_record, best_score

def extract_breed(description):
    match = re.search(r"\[breed=(.*?)\]", description)
    if match:
        return match.group(1)
    else:
        return "Unknown"


def main():

    args = parse_arguments()

    query_record = next(SeqIO.parse(args.query, "fasta"))
    database_records = list(SeqIO.parse(args.db, "fasta"))

    query_seq = str(query_record.seq)

    best_record, score = find_best_match(query_seq, database_records)

    percent_identity = score / len(query_seq) * 100
    breed = extract_breed(best_record.description)

    print("Closest breed:", breed)
    print("Best match ID:", best_record.id)
    print("Alignment score:", int(score))
    print("Percent identity:", f"{percent_identity:.2f}%")


if __name__ == "__main__":
    main() 