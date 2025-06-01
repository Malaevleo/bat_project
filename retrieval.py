import Bio
from Bio import Entrez
# in case it fails and wants to file a report
Entrez.email = "your@email.com"

with open("secretedfactorsids.txt") as f: # Place txt file with the REFSEQ ids of genes 
    refseq_ids = f.read().strip().replace(':', '\n').split()

with open("homo_secrets.fasta", "w") as output_fasta: # Output fasta file
    for refseq_id in refseq_ids:
        try:
            handle = Entrez.efetch(db="protein", id=refseq_id, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            output_fasta.write(fasta_data)
            print(f"Retrieved: {refseq_id}")
        except Exception as e:
            print(f"Error retrieving {refseq_id}: {str(e)}")