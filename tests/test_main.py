import sys
import os
import subprocess
import pytest

# Ensure main module can be imported from parent directory
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from main import (sequence_identity, extract_breed, find_best_match, 
                  load_fasta_records, validate_query_sequence, calculate_p_value)

# Test sequence_sequence identity with typical DNA comparisons
@pytest.mark.parametrize("seq1, seq2, expected_matches, expected_valid",
                         [
                             ("ATCG", "ATCG", 4, 4),
                             ("ATCG", "ATGC", 2, 4),
                             ("AAAA", "TTTT", 0, 4),
                         ],
                         )

def test_sequence_identity(seq1, seq2, expected_matches, expected_valid):
    matches, valid = sequence_identity(seq1, seq2)

    # Check correct match count and valid positions
    assert matches == expected_matches
    assert valid == expected_valid

# Test behaviour when no valid DNA bases are present
def test_sequence_identity_no_valid_bases():
    matches, valid = sequence_identity("XXXX", "YYYY")

    assert matches == 0 
    assert valid == 0

# Test breed extraction from FASTA description
def test_extract_breed():

    description = "dg|123| [breed=boxer] Canis lupus familiaris"
    breed = extract_breed(description)

    assert breed == "boxer"

# Test fallback when breed is missing
def test_extract_breed_missing():
    description = "no breed info"
    assert extract_breed(description) == "Unknown"

# Dummy class to simulate SeqRecord objects
class DummyRecord:
    def __init__(self, seq, id="test", description=""):
        self.seq = seq
        self.id = id
        self.description = description

# Test best match selection from database
def test_find_best_match():

    query = "ATCG"

    db = [DummyRecord("ATCG", "seq1"), 
          DummyRecord("ATGC", "seq2"), 
          DummyRecord("TTTT", "seq3")]

    best, matches, valid, probs = find_best_match(query, db)

    # Ensure correct sequence is selected
    assert best.id == "seq1"
    assert matches == 4 
    assert valid == 4

# Test that highest similarity is selected in comparisons
def test_find_best_match_selects_highest_similarity():
    query = "ATCG"

    db = [DummyRecord("ATCG", "perfect"),
          DummyRecord("ATCG", "almost"),]
    
    best, matches, valid, probs = find_best_match(query, db)

    assert best.id == "perfect"

# Test behaviour when multiple sequences tie
def test_find_best_match_tie():
    query = "ATCG"

    db = [DummyRecord("ATCG", "seq1"),
          DummyRecord("ATCG", "seq2")]
    
    best, matches, valid, probs = find_best_match(query, db)
    
    assert best.id in ["seq1", "seq2"]

# Test valid DNA sequence passes validation
def test_validate_query_sequence_valid():
    record = DummyRecord("ATCG")
    seq = validate_query_sequence(record)
    assert seq == "ATCG"

# Test invalid DNA sequence raises error
def test_validate_query_sequence_invalid():
    record = DummyRecord("XXXX")

    with pytest.raises(ValueError):
        validate_query_sequence(record)

# Test p-value calculation correctness
def test_calculate_p_value():
    p = calculate_p_value(matches=4, database_size=10)

    expected = (0.25 ** 4) * 10

    assert isinstance(p, float)
    assert p == expected

# Test sorting of similarity probabilities
def test_probabilities_sorted():
    probs = [("a", 0.5), ("b", 0.9), ("c", 0.7)]

    probs.sort(key=lambda x: x[1], reverse=True)

    assert probs[0][1] >= probs[1][1] >= probs[2][1]

# Test FASTA loading with valid file
def test_load_fasta_records(tmp_path):
    fasta_file = tmp_path / "test.fa"

    fasta_file.write_text(">seq1\nATCG\n")

    records = load_fasta_records(str(fasta_file))

    assert len(records) == 1
    assert str(records[0].seq) == "ATCG"

# Test FASTA loading with empty file raises error
def test_load_fasta_records_empty(tmp_path):
    fasta_file = tmp_path / "empty.fa"
    fasta_file.write_text("")

    with pytest.raises(ValueError):
        load_fasta_records(str(fasta_file))

# Integration test: program runs successfully
def test_program_runs():
    result = subprocess.run(
        ["python3", "main.py", "--db", "dog_breeds.fa", "--query", "mystery.fa"], 
        capture_output=True, text=True)

    assert result.returncode == 0
    assert "BEST MATCH" in result.stdout
    assert "SIMILARITY ACROSS ALL SEQUENCES" in result.stdout
    assert "P-value" in result.stdout

# Integration test: program runs with phylogeny option
def test_program_runs_with_phylogeny():
    result = subprocess.run(["python3", "main.py", "--db", "dog_breeds.fa", "--query", "mystery.fa", "--phylogeny"],
                            capture_output=True, text=True)
    
    assert result.returncode == 0
    assert "PHYLOGENETIC" in result.stdout
    



