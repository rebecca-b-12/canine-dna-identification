import sys
import os
import subprocess

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from main import sequence_identity, extract_breed, find_best_match

def test_sequence_identity_exact_match():

    seq1 = "ATCG"
    seq2 = "ATCG"

    matches, valid = sequence_identity(seq1, seq2)

    assert matches == 4
    assert valid == 4

def test_ssequence_identity_partial_match():

    seq1 = "ATCG"
    seq2 = "ATGC"

    matches, valid = sequence_identity(seq1, seq2)

    assert matches == 2
    assert valid == 4

def test_extract_breed():

    description = "dg|123| [breed=boxer] Canis lupus familiaris"
    breed = extract_breed(description)

    assert breed == "boxer"

class DummyRecord:
    def __init__(self, seq, id="test"):
        self.seq = seq
        self.id = id

def test_find_best_match():

    query = "ATCG"

    db = [DummyRecord("ATCG", "seq1"), DummyRecord("ATGC", "seq2"), DummyRecord("TTTT", "seq3")]

    best, matches, valid, probs = find_best_match(query, db)

    assert best.id == "seq1"

def test_program_runs():
    result = subprocess.run(["python", "main.py", "--db", "dog_breeds.fa", "--query", "mystery.fa"], capture_output=True, text=True)

    assert result.returncode == 0
    assert "Closest sequence" in result.stdout

    



