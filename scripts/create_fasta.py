import os
from Bio import SeqIO


def create_fasta(sequences, output, filename):
    """
    Procedure that receives a dictionary with the sequences of the intergenic regions and creates a FASTA file with the sequences.
    Prepends query id$ (dictionary key) to each sequence record ID
    
    Parameters:
        sequences (dict): a dictionary containing the sequences of the intergenic regions, with each query protein as key.
        output (str): the path to the output directory where MEME results will be saved.
    """
    
    #create list of sequences
    records = []
    for id, seqs in sequences.items():
        for seq in seqs:
            records.append(seq)
    
    output_file_path = os.path.join(output,filename)

    try:
        with open(output_file_path, mode="w") as output_fasta:
            SeqIO.write(records, output_fasta, "fasta")
    except Exception as e:
        print("Error while saving sequences:", e)    