from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from meme_aux import max_motif
import os
from logger import logger


def meme_motif_coverage(output, meme_output = False):
    ''' Function to compute the query coverage of a MEME motif
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
        
        The coverage is computed in a double loop, for MEME motif and for query within.
        The loop counts the number of instances associated with the motif for that particular query, avoiding double counts
        (i.e. when two instances are in the same sequence).
        It uses these counts to compute the micro and macro averages (over queries) for coverage, and also an integer coverage ratio.
        The micro-average is simply the total number of instances for the motif in sequences mapping to all queries,
        divided by the total number of sequences mapping to that motif.
        The macro-average is the average of the coverages per query (i.e. sum of coverages divided by # of queries).
        
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
        
    meme_coverage_path = os.path.join(output,"meme_coverage.csv")
    csv_file = open(meme_coverage_path,'w')

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
    csv_header1 = 'MotifN,Motif,Evalue,C_Fraction,C_microAvg,C_macroAvg'
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

    #dictionary with motif IDs as keys, and query coverage as values
    motif_stats = {}

    #compute query coverage for each motif
    #for each motif reported by MEME
    for meme_motif in meme_out:
        #define dictionary with keys set for all queries, 
        #the instance count dictionary intialized to zero, and
        #the total number of sequences per query set to what is read from input file
        motif_stats[meme_motif.name] = {}
        for query in list_of_queries:
            motif_stats[meme_motif.name][query] = {'instance_cnt' : 0, 'total_seqs' : queries[query], 'coverage' : None, 'int_coverage' : 0}

        #set of nucleotide IDs already processed
        #[used to check that multiple instances in same sequence are not double counted]
        used_seqs = set()
        #go through every motif instance
        for m_instance in meme_motif.instances:
            #obtain the query ID for the sequence in which the instance is
            m_query_id = m_instance.sequence_name.split('$')[0]
            #add to query tally
            #only add instance to count if the instance was NOT in a sequence already processed
            if m_instance.sequence_name not in used_seqs:
                motif_stats[meme_motif.name][m_query_id]['instance_cnt'] = motif_stats[meme_motif.name][m_query_id]['instance_cnt'] + 1
            #update last nucleotide id
            used_seqs.add(m_instance.sequence_name)
        #compute coverage for this motif on this query
        for query in list_of_queries:
            motif_stats[meme_motif.name][query]['coverage'] = motif_stats[meme_motif.name][query]['instance_cnt'] / motif_stats[meme_motif.name][query]['total_seqs']
            if motif_stats[meme_motif.name][query]['coverage'] > 0.0:
                motif_stats[meme_motif.name][query]['int_coverage'] = 1

        #add dictionary entries for micro and macro averages of coverage
        motif_stats[meme_motif.name]['coverage (int)']=0
        motif_stats[meme_motif.name]['coverage (uA)']=0
        motif_stats[meme_motif.name]['coverage (MA)']=0

        #compute micro and macro averages for coverage
        #micro average (sum of all instances over all total sequences, for all queries)
        #macro average (sum of all query coverages, over number of queries) 
        #int coverage (integer fraction)
        for query in list_of_queries:
            motif_stats[meme_motif.name]['coverage (int)'] = motif_stats[meme_motif.name]['coverage (int)'] + \
                                                             motif_stats[meme_motif.name][query]['int_coverage']
            motif_stats[meme_motif.name]['coverage (uA)'] =  motif_stats[meme_motif.name]['coverage (uA)'] + \
                                                             motif_stats[meme_motif.name][query]['coverage'] * motif_stats[meme_motif.name][query]['total_seqs'] / len(meme_input_seqs)
            motif_stats[meme_motif.name]['coverage (MA)'] =  motif_stats[meme_motif.name]['coverage (MA)'] + \
                                                             motif_stats[meme_motif.name][query]['coverage']
        motif_stats[meme_motif.name]['coverage (int)'] = motif_stats[meme_motif.name]['coverage (int)'] / len(list_of_queries)
        motif_stats[meme_motif.name]['coverage (MA)'] = motif_stats[meme_motif.name]['coverage (MA)'] / len(list_of_queries)

        #write to console
        #motif overall stats
        print('Motif:', meme_motif.name)
        print('|--> Coverage (int. fraction):', motif_stats[meme_motif.name]['coverage (int)'])
        print('|--> Coverage (micro average):', motif_stats[meme_motif.name]['coverage (uA)'])
        print('|--> Coverage (macro average):', motif_stats[meme_motif.name]['coverage (MA)'])
        logger.info('Motif: ' + meme_motif.name)
        logger.info('|--> Coverage (int. fraction): ' + str(motif_stats[meme_motif.name]['coverage (int)']))
        logger.info('|--> Coverage (micro average): ' + str(motif_stats[meme_motif.name]['coverage (uA)']))
        logger.info('|--> Coverage (macro average): ' + str(motif_stats[meme_motif.name]['coverage (MA)']))
        print('|--> Coverage, per query:')
        logger.info('|--> Coverage, per query:')        
        #query-specific motif stats
        for query in list_of_queries:
            print('    |--> Query:',query)
            print('         |--> Instances: ', motif_stats[meme_motif.name][query]['instance_cnt'])
            print('         |--> Sequences: ', motif_stats[meme_motif.name][query]['total_seqs'])
            print('         |--> Coverage:  ', motif_stats[meme_motif.name][query]['coverage'])
            logger.info('    |--> Query: ' + query)
            logger.info('         |--> Instances: ' + str(motif_stats[meme_motif.name][query]['instance_cnt']))
            logger.info('         |--> Sequences: ' + str(motif_stats[meme_motif.name][query]['total_seqs']))
            logger.info('         |--> Coverage:  ' + str(motif_stats[meme_motif.name][query]['coverage']))

        #write to CSV  
        #motif overall stats
        csv_text = meme_motif.alt_id + ',' + meme_motif.name + ',' + str(meme_motif.evalue) + ',' +\
                                           str(motif_stats[meme_motif.name]['coverage (int)']) + ',' +\
                                           str(motif_stats[meme_motif.name]['coverage (uA)']) + ',' +\
                                           str(motif_stats[meme_motif.name]['coverage (MA)'])
        #query-specific motif stats: coverage values
        for query in list_of_queries:
            csv_text = csv_text + ',' + str(motif_stats[meme_motif.name][query]['coverage'])
        #query-specific motif stats: ratio values
        for query in list_of_queries:
            csv_text = csv_text + ',' + str(motif_stats[meme_motif.name][query]['instance_cnt']) + '|' + \
                                        str(motif_stats[meme_motif.name][query]['total_seqs'])

        csv_file.write(csv_text + '\n')

    #write to console
    print('Best motif (int. fraction):', max_motif(motif_stats,'coverage (int)'))
    print('Best motif (micro average):', max_motif(motif_stats,'coverage (uA)'))
    print('Best motif (macro average):', max_motif(motif_stats,'coverage (MA)'))
    logger.info('Best motif (int. fraction): ' + str(max_motif(motif_stats,'coverage (int)')))
    logger.info('Best motif (micro average): ' + str(max_motif(motif_stats,'coverage (uA)')))
    logger.info('Best motif (macro average): ' + str(max_motif(motif_stats,'coverage (MA)')))
    
    csv_file.close()
    return(motif_stats)

