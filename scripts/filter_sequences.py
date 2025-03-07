from Bio import Align
from calculate_identity import calculate_identity
from logger import logger
    
    
def filter_sequences(sequences, configParams):
    """
    Filter sequences based on pairwise sequence identity.

    This function performs pairwise sequence alignment between the input sequences
    and filters out those that have a sequence identity greater than or equal to
    the specified maximum identity threshold.

    Parameters:
        sequences (list): a list of dictionaries containing sequence information.
            Each dictionary should contain a key "seq" with a BioPython SeqRecord object
            representing the sequence.
        configParams (Variables): an instance of the Variables class containing configuration parameters.
            The configuration parameters should include the maximum identity threshold.

    Returns:
        tuple: a tuple containing two elements:
            - A list of dictionaries containing the filtered sequences that meet the identity threshold.
            - A dictionary containing pairwise sequence identities for the retained sequences.
              The keys are sequence IDs, and the values are lists of tuples.
              Each tuple contains the ID of the aligned sequence, the identity score, and the aligned sequence itself.

    Note:
        The filtering process is performed by pairwise sequence alignment using the Needleman-Wunsch algorithm.
        Sequences with an identity greater than or equal to the maximum identity threshold are filtered out.
        The output retains only sequences that meet the identity threshold, along with their pairwise sequence identities.
        
        The algorithm is very simple. A double loop goes first through primary sequences and then through secondary ones.
        For each primary sequence, all secondary sequences are compared to it. Any showing higher than threshold identity
        are marked for deletion. The following primary sequence starts at the first sequence not marked for deletion,
        and omits comparisons with any sequences tagged for deletion.
    """
    
    max_identity = configParams.max_identity
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    idx_to_remove = set()
    identities_dict = {}
    
    # Flatten the nested list of sequences
    sequences_list = [d for subseq in sequences for d in subseq]
    
    # Check all sequences for validity, mark as to-delete if not valid
    idx_to_remove = set()
    for idx, sequence in enumerate(sequences_list):
         if len(sequence) == 0 or "N" in sequence:
            idx_to_remove.add(idx)
            logger.info(f"-Sequence {sequence.id} is not valid and has been removed from dataset")
    
    #report unfiltered sequences that are valid
    all_sequences = [seq for idx, seq in enumerate(sequences_list) if idx not in idx_to_remove]
    
    #initiate the 2-tier for loop for sequence filtering
    for idx, sequence in enumerate(sequences_list):
        # skip sequence if marked for deletion
        if idx in idx_to_remove:
            continue
            
        seq1 = str(sequence.seq)
        # If seq1 is empty or contains N, we mark it for deletion
        # seq1 = str(sequence.seq)
        # if len(seq1) == 0 or "N" in seq1:
            # idx_to_remove.add(idx)
            # continue 
        
        identities = []
        for idx2, sequence2 in enumerate(sequences_list[(idx + 1):]):
            # Get the real index of the sequence
            idx2 = idx2 + idx + 1
            # skip sequence if marked for deletion
            if idx2 in idx_to_remove:
                continue
            seq2 = str(sequence2.seq)
            # If seq2 is empty or contains N, we mark it for deletion
            # if len(seq2) == 0 or "N" in seq2:
                # idx_to_remove.add(idx2)
                # continue
            # If the two sequences are literally equal, mark seq2 for deletion
            if seq1 == seq2:
                idx_to_remove.add(idx2)
                continue

            # perform pairwise alignment between seq1 and seq2
            alignments = aligner.align(seq1, seq2)
            # compute the identity between both sequences
            identity = calculate_identity(alignments)

            # mark for deletion seq2 if identity above threshold
            if identity >= max_identity:
                idx_to_remove.add(idx2)
            # if identity is low enough, store it in identities dictionary
            # so that it can be used later to estimate divergence among motif instances
            else:
                identities.append((sequence.id, sequence2.id, identity))
        # once all secondary sequences have been processed, add identity list to dictionary for seq1
        identities_dict[sequence.id] = identities

    # Filter out the sequences to be discarded, leaving only sequences to keep
    sequences_to_keep = [seq for idx, seq in enumerate(sequences_list) if idx not in idx_to_remove]

    # print("Seqs to keep\n", sequences_to_keep)
    # print("All seqs\n", all_sequences)
    
    return all_sequences, sequences_to_keep, identities_dict
