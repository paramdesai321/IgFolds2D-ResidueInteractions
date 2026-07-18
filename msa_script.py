from Bio import SeqIO

# Provide the path to your FASTA file
fasta_file = "P01732_P10966_P01730_1RHH_MSA_merged.fasta"

# Read through the file and print sequence IDs and lengths
for record in SeqIO.parse(fasta_file, "fasta"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}")
    print(f"Length: {len(record.seq)}")

