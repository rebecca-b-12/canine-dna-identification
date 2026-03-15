import argparse
import re
from Bio import SeqIO
from Bio.Align import PairwiseAligner

def parse_arguments():
    parser = argparse.ArgumentParser(description="DNA Identification Service")
    parser.add_argument("--db", required=True, help="FASTA database file")
    parser.add_argument("--query", required=True, help="Query FASTA file")
    return parser.parse_args()

def find_best_match(query_seq, database_records):

    aligner = PairwiseAligner()
    aligner.mode = "global"

    best_record = None
    best_score = -1

    for record in database_records:

        score = aligner.score(query_seq, record.seq)

        if score > best_score:
            best_score = score
            best_record = record
    return best_record, best_score

def extract_breed(description):
    match = re.search(r"\[breed=(.*?)\]", description)
    if match:
        return match.group(1)
    return "Unknown"


def main():

    args = parse_arguments()

    query_record = next(SeqIO.parse(args.query, "fasta"))
    database_records = list(SeqIO.parse(args.db, "fasta"))

    best_record, score = find_best_match(query_record.seq, database_records)

    identity = score / len(query_record.seq) * 100

    breed = extract_breed(best_record.description)

    print("Closest breed:", breed)
    print("Best match ID:", best_record.id)
    print("Alignment score:", int(score))
    print("Percent identity:", f"{identity:.2f}%")


if __name__ == "__main__":
    main() 