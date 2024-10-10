# bat_project


Blastp results for Myotis Myotis (resultsmyofiltered.txt):
makeblastdb -in myotis_protein.faa -dbtype prot -out myotis_myotis_proteome
blastp -query matrisome/refseq_sequences.fasta -db myotis_myotis_proteome -out resultsmyofiltered.txt -outfmt 6 -num_threads 6 -max_target_seqs 3
