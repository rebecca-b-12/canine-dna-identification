# canine-dna-identification
A DNA-based identification service for assigning a query sequence to the most similar reference sequence in a canine genomic database, with statistical confidence estimation and phylogenetic reconstruction.

# Features
- Reads DNA sequences from FASTA files using Biopython
- Compares a query sequence to sequences in a database
- Identifies the closest matching sequence
- Computes similarity probabilities across the database
- Estimates a p-value for the match
- Includes unit tests and integration tests using pytest

# Requirements
- Python 3.10 or higher
- Biopython
- pytest (for running tests)
Install dependencies with:
pip install -r requirements.txt

Running the program
The program requires two inputs:
- A FASTA file containing known dog DNA sequences (database)
- A FASTA file containing an unkown query sequence
Example command:
python main.py --db dog_breeds.fa --query mystery.fa




# Example Output
Closest sequence: gb|AY656744.1|
Breed: English Springer Spaniel
Difference: 73

Probabilities across database:
Rank  Sequence ID           Probability
----------------------------------------
1     gb|AY656744.1|        0.9957
2     gb|CM023446.1|        0.9923
3     gb|AB123456.1|        0.9911

p_value: 0.0


# Running Tests
Unit and integration tests are included and can be run with: pytest
Tests verify:
- Sequence comparison logic
- Correct extraction of breed information
- Successful execution of the full program


# Project Structure
canine-dna-identification/
│
├── main.py                # main DNA identification program
├── dog_breeds.fa          # FASTA database of dog sequences
├── mystery.fa             # query sequence to identify
├── requirements.txt       # project dependencies
├── README.md              # project documentation
│
├── tests/
│   └── test_main.py       # unit and integration tests
│
└── venv/                  # virtual environment (optional)


# Author
Rebecca Bell
