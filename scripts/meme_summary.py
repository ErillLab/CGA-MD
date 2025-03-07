#from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs

import os

def meme_motif_summary(coverage_stats, divergence_stats, output, meme_output = False):
    '''Function that generates a summary CSV file with the coverage and divergence stats for a motif.

       Receives:
           - coverage_stats : the stats dictionary for coverage
           - divergence_stats : the stats dictionary for divergence
           - path to input file for MEME in FASTA format (the pipeline output folder)
           - path to MEME output file (the folder with results for the MEME analysis), if not within preespecified location

        Generates output CSV file.
    '''
    #set output folder for MEME if only rerunning MEME
    if meme_output: 
        meme_out_path = os.path.join(output,meme_output)
    else:
        meme_out_path = os.path.join(output,'meme_out')
        
    #open CSV for output
    meme_summary_path = os.path.join(output,"meme_summary.csv")
    csv_file = open(meme_summary_path,'w')
    
    #write header for CSV file
    csv_header = 'MotifN,Motif,Evalue,C_Fraction,C_microAvg,C_macroAvg,Site_uAvg,Site_MAvg,Seq_uAvg,Seq_MAvg\n'
    csv_file.write(csv_header)

    #load up meme output data
    try:
        with open(os.path.join(meme_out_path, 'meme.xml')) as meme_file:
            meme_out = motifs.parse(meme_file,'MEME')
    except Exception as e:
        print("Error reading MEME XML file:", e)
        return None

    #compute query-based divergence for each motif
    #for each motif reported by MEME
    for meme_motif in meme_out:
        #write to CSV  
        #motif overall stats
        csv_text = meme_motif.alt_id + ',' + meme_motif.name + ',' + str(meme_motif.evalue) + ',' +\
                                           str(coverage_stats[meme_motif.name]['coverage (int)']) + ',' +\
                                           str(coverage_stats[meme_motif.name]['coverage (uA)']) + ',' +\
                                           str(coverage_stats[meme_motif.name]['coverage (MA)']) + ',' +\
                                           str(divergence_stats[meme_motif.name]['site_divergence (uA)']) + ',' +\
                                           str(divergence_stats[meme_motif.name]['site_divergence (MA)']) + ',' +\
                                           str(divergence_stats[meme_motif.name]['seq_divergence (uA)']) + ',' +\
                                           str(divergence_stats[meme_motif.name]['seq_divergence (MA)'])

        csv_file.write(csv_text + '\n')
    csv_file.close()
    

