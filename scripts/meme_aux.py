'''Auxiliary functions for MEME analysis
'''

def max_motif(motif_stats,stat_type):
    '''Function to compute the best motif according to a statistic type
        - Inputs:
            - The motif_stats dictionary, with motif name as key and the statistic
              values as values (and a dictionary for each query statistic values for the motif)
            - The statistic type:
                - coverage (int)
                - coverage (uA)
                - coverage (MA)
                - 'site_divergence (uA)
                - 'site_divergence (MA)
                - 'seq_divergence (uA)'
                - 'seq_divergence (MA)'
                
        - Returns
            - Tuple with best motif name and its statitsic value
    '''
    max_stat = 0
    max_motif = ''
    for motif, data in motif_stats.items():
        if data[stat_type] > max_stat:
            max_stat = data[stat_type]
            max_motif = motif
    return (max_motif, max_stat)
    
    
def extend_site(meme_site,input_seqs, extension):
    '''Returns extended site with flanking regions, as well as the sequence index,
       using the site sequence id to locate the corresponding input sequence and extend on it
    '''
    #extend the instance to include left and right flanks (+/-10)
    #get site parameters
    seq_id = meme_site.sequence_name
    site_start = meme_site.start
    site_length = meme_site.length

    #get left and right flanks by accessing the original sequence
    #first identify corresponding sequence among input sequences
    seq_to_extend = None
    sequence_index = 0
    for idx, input_seq in enumerate(input_seqs):
        if (input_seq.id == seq_id):
            seq_to_extend = input_seq
            sequence_index = idx

    #obtain flanks
    left_flank = seq_to_extend[site_start-extension-1:site_start-1]
    right_flank = seq_to_extend[site_start+site_length-1:site_start+site_length+extension-1]
    center = seq_to_extend[site_start-1:site_start+site_length-1]
    site_seq = left_flank + center + right_flank
    #reverse-complement sequence if site was reported as reverse complement
    if meme_site.strand != '+':
        site_seq=site_seq.reverse_complement(id=site_seq.id,description=site_seq.description,name=site_seq.name)

    return(sequence_index, site_seq)
