def calculate_identity(alignment):
    """
    Calculates the percentage sequence identity of a sequence alignment.
    
    Parameters:
        alignment (Bio.Align.PairwiseAlignment): a Biopython pairwise alignment object,
            containing the alignment information between two sequences.

    Return:
        identity (float): a float value representing the percentage identity between the two aligned sequences.
            This value will be in the range of 0 to 1, where 1 indicates 100% identity.
    """
    
    symbol_line = format(alignment[0]).split("\n")[1]
    
    identical_count = symbol_line.count("|")
    
    identity = identical_count/len(symbol_line)
    
    return identity