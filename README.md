# canine-dna-identification
A DNA-based identification service for assigning a query sequence to the most similar reference sequence in a canine genomic database, with statistical confidence estimation and phylogenetic reconstruction.

What main.py should do
1. parse command line arguments
2. load sequences from a FASTA database
3. Take a query sequence
4. Align query to every database sequence
5. Find best match
6. Print: best sequence ID, alignment score and percent identity