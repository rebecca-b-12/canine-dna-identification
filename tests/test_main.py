import sys
import os
import subprocess
import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from main import sequence_identity, extract_breed, find_best_match

@pytest.mark.parametrize("seq1, seq2, expected_matches, expected_valid",
                         [("ATCG", "ATCG", 4, 4),
                          ("ATCG", "ATGC", 2, 4),
                          ("AAAA", "TTTT", 0, 4),],)

def test_sequence_identity(seq1, seq2, expected_matches, expected_valid):
    matches, valid = sequence_identity(seq1, seq2)

    assert matches == expected_matches
    assert valid == expected_valid

def test_extract_breed():

    description = "dg|123| [breed=boxer] Canis lupus familiaris"
    breed = extract_breed(description)

    assert breed == "boxer"

class DummyRecord:
    def __init__(self, seq, id="test", description=""):
        self.seq = seq
        self.id = id
        self.description = description

def test_find_best_match():

    query = "ATCG"

    db = [DummyRecord("ATCG", "seq1"), DummyRecord("ATGC", "seq2"), DummyRecord("TTTT", "seq3")]

    best, matches, valid, probs = find_best_match(query, db)

    assert best.id == "seq1"

def test_best_match_highest_similarity_selected():
    query = "ATCG"

    db = [DummyRecord("ATCG", "perfect"),
          DummyRecord("ATCG", "almost"),]
    
    best, matches, valid, probs = find_best_match(query, db)

    assert best.id == "perfect"

def test_program_runs():
    result = subprocess.run(["python3", "main.py", "--db", "dog_breeds.fa", "--query", "mystery.fa"], capture_output=True, text=True)

    assert result.returncode == 0
    assert "Closest sequence" in result.stdout
    assert "Similarity to query sequence" in result.stdout

    



