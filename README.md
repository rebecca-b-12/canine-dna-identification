# canine-dna-identification
A DNA-based identification tool for assigning a query sequence to the most similar reference sequence in a canine genomic database, with statistical confidence estimation and optional phylogenetic reconstruction.

# Features
- Reads DNA sequences from FASTA files using Biopython
- Compares a query sequence against a database of sequences
- Identifies the closest matching sequence
- Computes similarity scores across sequences
- Calculates the number of differing bases
- Estimates statistical significance using a p-value
- Displays results in a clean, formatted table
- Generates phylogenetic trees
- Highlights the best match in outputs
- Includes unit tests and integration tests using pytest

# Requirements
- Python 3.10 or higher
- Biopython
- pandas
- matplotlib
- tabulate
- pytest (for running tests)

Install dependencies with:
pip install -r requirements.txt

# Running the program
The program requires two inputs:
- A FASTA file containing known dog DNA sequences (database)
- A FASTA file containing an unknown query sequence
Basic usage:
- python main.py --db data/dog_breeds.fa --query data/mystery.fa
With phylogenetic analysis:
- python main.py --db data/dog_breeds.fa --query data/mystery.fa --phylogeny


# Example Output
Best Match

    Closest sequence: gb|AY656744.1|
    Breed: English Springer Spaniel
    Difference (number of differing bases): 73

Similarity Table
- Displays all sequences ranked by similarity
- Includes p-values for statistical confidence
- Best match is clearly highlighted

Summary
    Best match similarity: 0.9994
    Best match p-value: <1e-300

# Phylogenetic Analysis
When the --phylogeny flag is used:
- A full phylogenetic tree (including all sequences) is printed in the terminal
- A cleaner tree image is saved as: phylogenetic_tree.png
 
 Notes:
 - Method: Neighbour-Joining
 - Distance metric: 1 - sequence similarity
 - The image includes the top most similar sequences for readability
 - The full dataset is used for analysis

# Explanation of Output
- Similarity: proportion of matching DNA bases (range: 0-1)
- Difference: number of positions where the sequences differ
- P-value: estimated probability that the observed similarity occurred by chance

# Running Tests
Unit and integration tests are included and can be run with: pytest

Tests verify:
- Sequence comparison logic
- Best match identification
- Breed extraction
- P-value calculation
- Full program execution

# Project Structure
canine-dna-identification/
│
├── main.py              # main DNA identification program
├── requirements.txt     # project dependencies
├── README.md            # project documentation
│
├── data/
│   ├── dog_breeds.fa    # FASTA database of dog sequences
│   └── mystery.fa       # query sequence
│
├── tests/
│   └── test_main.py     # unit and integration tests
│
└── venv/                # virtual environment (optional)

# Assumptions
- DNA bases are assumed to have equal probability (0.25 each)
- Only valid bases (A, T, C, G) are considered in comparisons
- P-values are estimated using a simplified random model

# Future Improvements
- Use more advanced statistical models (e.g. binomial test)
- Improve phylogenetic tree visualisation
- Support multiple query sequences
- Add graphical user interface (GUI)

# Author
Rebecca Bell
