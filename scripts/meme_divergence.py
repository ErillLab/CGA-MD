from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio import Align
from math import isnan
from calculate_identity import calculate_identity
from meme_aux import max_motif
from meme_aux import extend_site
from logger import logger

import os


def meme_motif_divergence(output, meme_output = False):
    ''' Function to compute the site and sequence divergence of a MEME motif
        - Inputs:
            - path to input file for MEME in FASTA format (the pipeline output folder)
            - path to MEME output file (the folder with results for the MEME analysis), if not within preespecified location

        Headers in the input file for MEME is assumed to have the following format:
            Query_ID $ Ortho_ID1 | ... | OrthoIDN @ Nucleotide_ID
            where query_ID is the query ID, separated from orthologs by $
            the ortholog protein IDs are separated by pipes
            and after @ follows the nucleotide record with the genes encoding those orthologs

            Example:
            >WP_000154162.1$WP_053038396.1@NZ_CTXC01000075.1 Region between WP_000154162.1$WP_053038396.1 in NZ_CTXC01000075.1
            AACAATCATTTCCTTTCAAACAACGCCTTCTTTATATTTTAATTAGTGTCACTATTATTC
            ATTGTTTAAATTCCCTCATCTTAACTTATAGAAACATTTCCTCAGTGCAAATACTGTGAT
            GCTCAAATTAGTTGAACTGCTTTCCAATAATCCAAATTCACCCGAGAAAATCACGATATA
            TTAATATCTTTCAACTATGTAAGATTCAAACATTGAATTGTAAGTTTAGATAAAGGCTAC
            TGTTCTACTTATTATTTAAAATGAATTTAAATACAAAAGGAGTGATACT
            
        The function first loads the MEME input (FASTA) file and the MEME output (XML) file, which is parsed
        with the Bio.Motifs parser.
        
        From the FASTA input file, it gets all the query IDs and the number of ortholog sequences associated with them.
        
        The divergence computation is computed in a double loop, for MEME motif and for query within.
        Within this loop, another double loop evaluates the divergence among all sequence pairs mapping to instances of 
        the motif in sequences mapping to query.

        The evaluated divergences are used to compute the micro and macro averages (over queries) for site and sequence divergence.
        The micro-average is simply the sum of divergences for the motif in sequences mapping to all queries, divided by
        the total number of pairwise comparisons performed.
        The macro-average is the average of the divergences per query (i.e. sum of divergences divided by # of queries).
        
        All these values are stored in a dictionary that has the motif name as key, the average statistics as values,
        and a value that is a dictionary with the query IDs as keys and the specific statistics for that motif-query pair
        as values.
    '''
    #define location of output FASTA file
    meme_in_path = os.path.join(output,"sequences_filtered.fasta")

    #set output folder for MEME if only rerunning MEME
    if meme_output: 
        meme_out_path = os.path.join(output,meme_output)
    else:
        meme_out_path = os.path.join(output,'meme_out')

    meme_diversity_path = os.path.join(output,"meme_diversity.csv")
    csv_file = open(meme_diversity_path,'w')

    #dictionary with query IDs as key and number of ortholog sequences as value
    queries = {}

    #load up meme input sequence data
    #input file is assumed to be:
    # FASTA format
    # header: QUERY_ID$ortholog_id_list@nucleotide_record_id
    meme_input_seqs = list(SeqIO.parse(meme_in_path, 'fasta'))

    #obtain all query IDs and the number of ortholog sequences for each
    #for each sequence in the input set for MEME
    for input_seq in meme_input_seqs:
        #obtain query ID from header
        query_ID = input_seq.id.split('$')[0] 
        #if we have processed the query before, add current sequence to count
        if query_ID in queries.keys():
            queries[query_ID] = queries[query_ID] + 1
        #otherwise, initialize count for query
        else:
            queries[query_ID] = 1

    list_of_queries = queries.keys()

    #write header for CSV file
    csv_header1 = 'MotifN,Motif,Evalue,Site_uAvg,Site_MAvg,Seq_uAvg,Seq_MAvg'
    csv_header2 = ''
    for query in list_of_queries:
        csv_header2 = csv_header2 + ', ' + query 
    csv_header = csv_header1 + csv_header2 + csv_header2 + '\n'
    csv_file.write(csv_header)

    #load up meme output data
    try:
        with open(os.path.join(meme_out_path, 'meme.xml')) as meme_file:
            meme_out = motifs.parse(meme_file,'MEME')
    except Exception as e:
        print("Error reading MEME XML file:", e)
        return None     

    #initalize sequence aligner object
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    #dictionary with motif IDs as keys, and query-based divergences as values
    motif_divergences = {}

    #compute query-based divergence for each motif
    #for each motif reported by MEME
    for meme_motif in meme_out:
        #define dictionary with keys set for all queries, 
        #the distance dictionary intialized to zero, and
        #the total number of sequences per query set to what is read from input file    
        motif_divergences[meme_motif.name] = {}
        for query in list_of_queries:
            motif_divergences[meme_motif.name][query] = {'site_distance_sum' : 0, 'seq_distance_sum' : 0, 'total_seqs' : queries[query], 'total_comparisons' : 0, 'site_divergence' : None, 'seq_divergence' : None}

        total_motif_instances = len(meme_motif.instances)
        #for each query
        for query in list_of_queries:
            #for each instance
            for i in range(total_motif_instances):
                #if the site instance maps to current query
                if meme_motif.instances[i].sequence_name.split('$')[0] == query:
                    #extend site with flanking regions
                    seq_index_i, i_seq=extend_site(meme_motif.instances[i],meme_input_seqs, 10)
                    for j in range(i + 1, total_motif_instances):
                        #if the second site instance maps to current query
                        if meme_motif.instances[j].sequence_name.split('$')[0] == query:                
                            #extend site with flanking regions
                            seq_index_j, j_seq=extend_site(meme_motif.instances[j],meme_input_seqs, 10)

                            #align the two extended sites and compute distance
                            site_alignment = aligner.align(i_seq.seq, j_seq.seq)
                            site_distance = 1 - calculate_identity(site_alignment)

                            #align the two entire sequences harboring the sites and compute distance
                            site_alignment = aligner.align(meme_input_seqs[seq_index_i].seq, meme_input_seqs[seq_index_j].seq)
                            seq_distance = 1 - calculate_identity(site_alignment)

                            #store into divergence and update # of comparisons
                            motif_divergences[meme_motif.name][query]['site_distance_sum'] = \
                                              motif_divergences[meme_motif.name][query]['site_distance_sum'] + site_distance
                            motif_divergences[meme_motif.name][query]['seq_distance_sum'] = \
                                              motif_divergences[meme_motif.name][query]['seq_distance_sum'] + seq_distance
                            motif_divergences[meme_motif.name][query]['total_comparisons'] = \
                                              motif_divergences[meme_motif.name][query]['total_comparisons'] + 1

            #compute site divergence for the query do so only if there is at least one comparison
            if motif_divergences[meme_motif.name][query]['total_comparisons'] > 0:
                motif_divergences[meme_motif.name][query]['site_divergence'] = \
                                  motif_divergences[meme_motif.name][query]['site_distance_sum'] / \
                                  motif_divergences[meme_motif.name][query]['total_comparisons']
            #otherwise, set as NaN
            else:
                motif_divergences[meme_motif.name][query]['site_divergence'] = float("nan")
            #compute sequence divergence for the query do so only if there is at least one comparison
            if motif_divergences[meme_motif.name][query]['total_comparisons'] > 0:
                motif_divergences[meme_motif.name][query]['seq_divergence'] = \
                                  motif_divergences[meme_motif.name][query]['seq_distance_sum'] / \
                                  motif_divergences[meme_motif.name][query]['total_comparisons']
            #otherwise, set as NaN
            else:
                motif_divergences[meme_motif.name][query]['seq_divergence'] = float("nan")                                                                       

        #add dictionary entries for micro and macro averages of coverage
        motif_divergences[meme_motif.name]['site_divergence (uA)']=0
        motif_divergences[meme_motif.name]['site_divergence (MA)']=0
        motif_divergences[meme_motif.name]['seq_divergence (uA)']=0
        motif_divergences[meme_motif.name]['seq_divergence (MA)']=0

        #compute micro and macro averages for coverage
        #micro average (sum of all instances over all total sequences, for all queries)
        #macro average (sum of all query coverages, over number of queries) 
        #int coverage (integer fraction)
        query_cnt = 0
        comp_cnt = 0
        for query in list_of_queries:
            #do not consider motif-query pairs where there is no distance (none or 1 instances)
            if not isnan(motif_divergences[meme_motif.name][query]['site_divergence']):
                #update site divergence counts for average
                motif_divergences[meme_motif.name]['site_divergence (uA)'] =  motif_divergences[meme_motif.name]['site_divergence (uA)'] + \
                                                                 motif_divergences[meme_motif.name][query]['site_distance_sum']
                motif_divergences[meme_motif.name]['site_divergence (MA)'] =  motif_divergences[meme_motif.name]['site_divergence (MA)'] + \
                                                                 motif_divergences[meme_motif.name][query]['site_divergence']
                #update sequence divergence counts for average
                motif_divergences[meme_motif.name]['seq_divergence (uA)'] =  motif_divergences[meme_motif.name]['seq_divergence (uA)'] + \
                                                                 motif_divergences[meme_motif.name][query]['seq_distance_sum']
                motif_divergences[meme_motif.name]['seq_divergence (MA)'] =  motif_divergences[meme_motif.name]['seq_divergence (MA)'] + \
                                                                 motif_divergences[meme_motif.name][query]['seq_divergence']
                query_cnt = query_cnt + 1
                comp_cnt = comp_cnt + motif_divergences[meme_motif.name][query]['total_comparisons']

        if query_cnt > 0:
            motif_divergences[meme_motif.name]['site_divergence (MA)'] = motif_divergences[meme_motif.name]['site_divergence (MA)'] / query_cnt
            motif_divergences[meme_motif.name]['site_divergence (uA)'] = motif_divergences[meme_motif.name]['site_divergence (uA)'] / comp_cnt
            motif_divergences[meme_motif.name]['seq_divergence (MA)'] = motif_divergences[meme_motif.name]['seq_divergence (MA)'] / query_cnt
            motif_divergences[meme_motif.name]['seq_divergence (uA)'] = motif_divergences[meme_motif.name]['seq_divergence (uA)'] / comp_cnt
        else:
            motif_divergences[meme_motif.name]['site_divergence (MA)'] = float('nan')
            motif_divergences[meme_motif.name]['site_divergence (uA)'] = float('nan')
            motif_divergences[meme_motif.name]['seq_divergence (MA)'] = float('nan')
            motif_divergences[meme_motif.name]['seq_divergence (uA)'] = float('nan')

        #write to console
        #motif overall stats
        print('Motif:', meme_motif.name)
        print('|--> Site divergence (micro average):     ', motif_divergences[meme_motif.name]['site_divergence (uA)'])
        print('|--> Site divergence  (macro average):    ', motif_divergences[meme_motif.name]['site_divergence (MA)'])
        print('|--> Sequence divergence (micro average): ', motif_divergences[meme_motif.name]['seq_divergence (uA)'])
        print('|--> Sequence divergence  (macro average):', motif_divergences[meme_motif.name]['seq_divergence (MA)'])
        logger.info('Motif: ' + meme_motif.name)
        logger.info('|--> Site divergence (micro average):     ' + str(motif_divergences[meme_motif.name]['site_divergence (uA)']))
        logger.info('|--> Site divergence  (macro average):    ' + str(motif_divergences[meme_motif.name]['site_divergence (MA)']))
        logger.info('|--> Sequence divergence (micro average): ' + str(motif_divergences[meme_motif.name]['seq_divergence (uA)']))
        logger.info('|--> Sequence divergence  (macro average):' + str(motif_divergences[meme_motif.name]['seq_divergence (MA)']))
        #query-specific motif stats
        for query in list_of_queries:
            print('    |--> Query:',query)
            print('         |--> Total seqs:  ', motif_divergences[meme_motif.name][query]['total_seqs'])
            print('         |--> Site diver:  ', motif_divergences[meme_motif.name][query]['site_divergence'])
            print('         |--> Seq diver:   ', motif_divergences[meme_motif.name][query]['seq_divergence'])
            print('         |--> Total comps: ', motif_divergences[meme_motif.name][query]['total_comparisons'])
            logger.info('    |--> Query: ' + query)
            logger.info('         |--> Total seqs:  ' + str(motif_divergences[meme_motif.name][query]['total_seqs']))
            logger.info('         |--> Site diver:  ' + str(motif_divergences[meme_motif.name][query]['site_divergence']))
            logger.info('         |--> Seq diver:   ' + str(motif_divergences[meme_motif.name][query]['seq_divergence']))
            logger.info('         |--> Total comps: ' + str(motif_divergences[meme_motif.name][query]['total_comparisons']))
            
        #write to CSV  
        #motif overall stats
        csv_text = meme_motif.alt_id + ',' + meme_motif.name + ',' + str(meme_motif.evalue) + ',' +\
                                           str(motif_divergences[meme_motif.name]['site_divergence (uA)']) + ',' +\
                                           str(motif_divergences[meme_motif.name]['site_divergence (MA)']) + ',' +\
                                           str(motif_divergences[meme_motif.name]['seq_divergence (uA)']) + ',' +\
                                           str(motif_divergences[meme_motif.name]['seq_divergence (MA)'])

        #query-specific motif stats: site values
        for query in list_of_queries:
            csv_text = csv_text + ',' + str(motif_divergences[meme_motif.name][query]['site_divergence'])
        #query-specific motif stats: seq values
        for query in list_of_queries:
            csv_text = csv_text + ',' + str(motif_divergences[meme_motif.name][query]['seq_divergence'])

        csv_file.write(csv_text + '\n')

    #write to console
    print('Best motif (site micro average):', max_motif(motif_divergences,'site_divergence (uA)'))
    print('Best motif (site macro average):', max_motif(motif_divergences,'site_divergence (MA)'))
    print('Best motif (seq micro average): ', max_motif(motif_divergences,'seq_divergence (uA)'))
    print('Best motif (seq macro average): ', max_motif(motif_divergences,'seq_divergence (MA)'))
    logger.info('Best motif (site micro average): ' + str(max_motif(motif_divergences,'site_divergence (uA)')))
    logger.info('Best motif (site macro average): ' + str(max_motif(motif_divergences,'site_divergence (MA)')))
    logger.info('Best motif (seq micro average):  ' + str(max_motif(motif_divergences,'seq_divergence (uA)')))
    logger.info('Best motif (seq macro average):  ' + str(max_motif(motif_divergences,'seq_divergence (MA)')))

    csv_file.close()
    return(motif_divergences)

