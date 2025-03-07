import os
import argparse
from tqdm import tqdm
from Bio import Entrez
from variables import *
from logger import logger
from search_blast import search_blast
from get_gene_record import get_gene_record
from filter_sequences import filter_sequences
from read_and_validate_config import read_config
from extract_params_config import extract_params
from launch_meme_search import launch_meme_search

from meme_coverage import meme_motif_coverage
from meme_divergence import meme_motif_divergence
from meme_summary import meme_motif_summary

from grab_intergenic import grab_intergenic_regions
from create_fasta import create_fasta

import pickle

def main(root):
    pathConfig = os.path.join(root,"config/config.json")
    
    config = read_config(pathConfig)

    configParams = extract_params(config)
    
    Entrez.email = configParams.EMAIL
    Entrez.api_key = configParams.API_KEY
    run_meme_only = configParams.run_only_meme

    print("    |--> Crawl-back: ", "activated" if configParams.crawl_back else "not activated")
    logger.info("    |--> Crawl-back: " + "activated" if configParams.crawl_back else "not activated")
    print("    |--> Fuse intergenic: ", "activated" if configParams.fuse_intergenic else "not activated")
    logger.info("    |--> Fuse intergenic: " + "activated" if configParams.fuse_intergenic else "not activated")
    
    output = os.path.join(root, config["output_parameters"]["folder_name"])
    
    if not os.path.isdir(output):
        os.mkdir(output)
        
    #initialize main variables for sequence collection
    unfiltered_sequences, filtered_sequences, total_sequences, identities_id_dict = {}, {}, {}, {}
    
    # running in normal mode (BLAST + sequence retrieval)
    if run_meme_only == False:
        meme_output = False
        
        #BLAST mode (run BLAST to obtain orthologs and process them)
        if configParams.max_hits:
            # initiate BLAST search. will return a dictionary with input record ID as key and a list of hit accession numbers as the value
            # a BLAST search is performed for every protein ID in the config file
            hits = search_blast(params=configParams)
            hitsfile = open('hits', 'wb')
            pickle.dump(hits, hitsfile)
            
            # hitsfile = open('hits', 'rb')
            # hits = pickle.load(hitsfile)
        #non-BLAST mode (process only the input protein records)
        else:
            print('Skipping BLAST search')
            logger.info('Skipping BLAST search')
            hits = {}    #fake hits dictionary, will contain input protein IDs as key and value 
            for protID in configParams.input_records:
                hits[protID] = [protID]
        
        # process each BLAST hit
        for query_id, hits in hits.items():
            #list that will store retrieved intergenic sequences for all hits of this query ID
            sequences = []
            print(f"|--> Obtaining the gene records for the {query_id} hits. See the log file for more information.")
            
            for hit in tqdm(hits):
                # obtain from NCBI (IPG) the best gene record associated with the orthologous protein
                # the gene record will be prioritized based on completeness of the genome
                # the returned record is a dictionary with the nucleotide accession, strand, start and end positions
                record = get_gene_record(hit, configParams)
                if record is None:
                    logger.info(f"Skipping hit {hit} due to failed record retrieval.")
                    continue
                # grab relevant intergenic regions upstream of ortholog-encoding gene
                # the function returns a list of SeqRecords belonging to relevant intergenic regions upstream of the
                # gene of interest (record), with each intergenic region annotated with the protein_id and locus_tag 
                # of the gene it lies upstream of
                sequences.append(grab_intergenic_regions(query_id,record,configParams))
                
            print("    |--> Filtering the sequences...")
            logger.info("    |--> Filtering the sequences...")
            # filters obtained intergenic regions (for THIS hit) according to sequence identity
            # only intergenic regions that are different enough are kept
            # this returns also a dictionary of identities to be used in estimating motif instance diversity
            sequences_unfiltered, sequences_filtered, identities_dict = filter_sequences(sequences, configParams)

            #store retrieved & filtered sequences [intergenic regions] for this query in dictionary using query ID as key
            filtered_sequences[query_id] = sequences_filtered
            unfiltered_sequences[query_id] = sequences_unfiltered
            
            print(f"    |--> Intergenic regions found for {query_id}: {len(unfiltered_sequences[query_id])}")
            logger.info(f"    |--> Intergenic regions found for {query_id}: {len(unfiltered_sequences[query_id])}")
            print(f"    |--> Intergenic regions after filtering for {query_id}: {len(filtered_sequences[query_id])}")
            logger.info(f"    |--> Intergenic regions after filtering for {query_id}: {len(filtered_sequences[query_id])}")            
            
            #store identities 
            identities_id_dict[query_id] = identities_dict
            total_sequences[query_id] = len(sequences_filtered)

        print("    |--> Saving sequences in FASTA format...")
        logger.info("    |--> Saving sequences in FASTA format...")
        #generate FASTA file with identified & filtered sequences for all hits
        create_fasta(filtered_sequences, output, "sequences_filtered.fasta")
        create_fasta(unfiltered_sequences, output, "sequences_unfiltered.fasta")
    
    #running only to repeat MEME on data already compiled
    else:
        meme_output = config["output_parameters"]["folder_rerun_meme"]
    
    # launch MEME motif discovery
    print("Running MEME... (this can take a while)")
    logger.info("Running MEME...")
    launch_meme_search(configParams, output, meme_output)
    logger.info("|--> MEME search finished.")
    print("|--> MEME search finished.")
    
    logger.info("Reading MEME output and creating output files...")
    print("Reading MEME output and creating output files...")

    #process MEME output files
    #compute coverage for MEME motifs
    logger.info("|--> Computing MEME motif coverage.")
    print("|--> Computing MEME motif coverage.")
    coverage_stats = meme_motif_coverage(output, meme_output)
    #compute divergence for MEME motifs
    logger.info("|--> Computing MEME motif divergence.")
    print("|--> Computing MEME motif divergence.")
    divergence_stats = meme_motif_divergence(output, meme_output)
    #write out summary
    logger.info("|--> Writing out MEME analysis summary.")
    print("|--> Writing out MEME analysis summary.")
    meme_motif_summary(coverage_stats, divergence_stats, output, meme_output)
    logger.info("Output files created.")
    print("Output files created.")
    
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the main pipeline with the specified root directory.")
    parser.add_argument("root", type=str, help="The root directory for the pipeline.")
    args = parser.parse_args()
    
    main(args.root)