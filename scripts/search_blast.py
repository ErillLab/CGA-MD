from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from tqdm import tqdm
from logger import logger
import time

def build_taxon_query(taxIDs):
    """
    Builds a taxon query string for searching in biological databases.

    Parameters:
        taxIDs (list of int): s list of taxonomy IDs (integers).

    Returns:
        str: a query string formatted for use in biological database searches.

    Description:
        This function takes a list of taxonomy IDs and constructs a query string 
        that can be used to search for organisms in biological databases like NCBI. 
        The format of the query string depends on the number of taxonomy IDs provided:

    - If there is only one taxonomy ID, the function returns a query string in the format:
        "txid<taxID>[ORGN]".

    - If there are multiple taxonomy IDs, the function returns a query string 
        where each taxonomy ID is formatted as "txid<taxID>[ORGN]" and separated 
        by " OR ".
    """
    
    if len(taxIDs) == 1:
        return "txid" + ",".join(map(str, taxIDs)) + "[ORGN]"
    else:
        return " OR ".join([f"txid{taxID}[ORGN]" for taxID in taxIDs])

def search_blast(params, db = "nr", min_coverage = 0.75, max_hits = 50, e_value = 10E-10, taxIDs = None):
    """
    Perform a BLAST search for a set of protein records.

    This function conducts a BLAST search for a set of protein records using the specified parameters.
    It retrieves protein sequences from NCBI using Entrez, performs BLAST searches via NCBIWWW,
    and filters the search results based on specified criteria.

    Parameters:
        params (Variables): an instance of the Variables class containing configuration parameters.
        db (str, optional): the name of the BLAST database to be queried. Defaults to "nr".
        max_hits (int, optional): the maximum number of hits to return. Defaults to 50.
        e_value (float, optional): the threshold for the E-value of the hits. Defaults to 10E-10.
        taxIDs (list of strings, optional): the taxonomic IDs limits to the BLAST search. Defaults to None.
        min_coverage (float, optional): the minimum coverage for the hits. Defaults to None.

    Returns:
        hits_dict (dict): A dictionary containing the input record ID as the key and a list of hit accession numbers as the value.
    
    Raises:
        Exception: if no protein ID is provided for the BLAST search.
        Exception: if the BLAST search fails after the maximum number of attempts (REQUEST_LIMIT).

    Note:
        This function logs various steps of the BLAST search process using the provided logger.
    """
    
    REQUEST_LIMIT = params.REQUEST_LIMIT
    SLEEP_TIME = params.SLEEP_TIME
    
    # protein id's to BLAST
    input = params.input_records
    
    # config parameters for BLAST
    db = params.db
    taxIDs = params.tax_IDs
    min_coverage = params.query_coverage
    max_hits = params.max_hits
    e_value = params.e_value
    
    print("Starting BLAST search:")
    logger.info("Starting BLAST search:")
    
    if len(input) < 1:
        raise Exception("|--> At least one protein ID is needed to perform BLAST. Revise config file and try again.")
    
    hits_dict = {}
    blast_records = None
    
    # for each protein accession in config file...
    for id in tqdm(input):
        print("|--> Getting protein record")
        logger.info("|--> Getting protein record")
        
        #  obtain the protein sequence
        for i in range(REQUEST_LIMIT):
            
            try:
                handle = Entrez.efetch("protein", id=id, rettype="fasta", retmode="text")
                time.sleep(SLEEP_TIME)
                break
                
            except:
                logger.warning(f"    |--> NCBI exception raised on attempt {str(i)}. Reattempting.")
                
                if i == (REQUEST_LIMIT - 1):
                    logger.info(f"    |--> Could not download record after {str(REQUEST_LIMIT)} attempts.")
                    raise Exception(f"    |--> Could not download record after {str(REQUEST_LIMIT)} attempts.")
        
        logger.info("|--> Getting protein sequence")
        input_seq = (SeqIO.read(handle, "fasta")).seq
        
        print(f"|--> Performing BLAST search: {str(id)}")
        logger.info(f"|--> Performing BLAST search: {str(id)}")
        
        # perform the BLAST search with the obtained protein sequence
        if taxIDs:
            # search restricted by taxon IDs from config file
            taxon_query = build_taxon_query(taxIDs)
            print(f"    |--> Performing BLAST search for {str(id)} on taxID {",".join(map(str, taxIDs))}")
            logger.info(f"    |--> Performing BLAST search on taxID {",".join(map(str, taxIDs))}")
            for i in range(REQUEST_LIMIT):
                try:
                    result_handle = NCBIWWW.qblast(program="blastp",
                                                   database=db, sequence=input_seq,
                                                   entrez_query=taxon_query, expect=e_value,
                                                   hitlist_size=max_hits)
                    logger.info("|--> Getting records.")
                    blast_records = list(NCBIXML.parse(result_handle))
                    time.sleep(SLEEP_TIME)
                    break
                    
                except:
                    print(f"    |--> NCBI exception raised on attempt {str(i)}. Reattempting.")
                    logger.warning(f"    |--> NCBI exception raised on attempt {str(i)}. Reattempting.")

                    if i == (REQUEST_LIMIT - 1):
                        print(f"    |--> Could not download record after {str(REQUEST_LIMIT)} attempts.")
                        logger.error(f"    |--> Could not download record after {str(REQUEST_LIMIT)} attempts.")
        else:
            # search not restricted by taxon IDs
            for i in range(REQUEST_LIMIT):
                logger.info(f"    |--> Performing unrestricted BLAST search for {str(id)}")
                try:
                    result_handle = NCBIWWW.qblast("blastp", db, input_seq, expect=e_value, hitlist_size=max_hits)
                    logger.info("|--> Getting records.")
                    blast_records = list(NCBIXML.parse(result_handle))
                    time.sleep(SLEEP_TIME)
                    break
                    
                except:
                    print(f"    |--> NCBI exception raised on attempt {str(i)}. Reattempting.")
                    logger.warning(f"    |--> NCBI exception raised on attempt {str(i)}. Reattempting.")
                
                if i == (REQUEST_LIMIT - 1):
                    print(f"    |--> Could not download record after {str(REQUEST_LIMIT)} attempts.")
                    logger.error(f"    |--> Could not download record after {str(REQUEST_LIMIT)} attempts.")
        
        hits_dict[id] = []
        
        # iterate through the obtained BLAST alignments
        for record in blast_records[0].alignments:
            
            # get the hit id
            current_hit = record.hit_id.split("|")[1]
            logger.info(f"    |--> Analyzing hit {str(current_hit)}")
            
            # check that hit meets required min coverage criterion from config file and store
            for hit in record.hsps:
                # analyze the coverage as the fraction of query aligned
                if min_coverage:
                    cov = (hit.query_end - hit.query_start + 1) / len(input_seq)
                    # if coverage is sound, append hit to hit list under dictionary entry for this query id
                    if cov >= min_coverage:
                        if current_hit not in hits_dict[id]:
                            logger.info(f"    |--> Adding hit (coverage = {str(round(cov*100, 2))}): {str(current_hit)}.")
                            hits_dict[id].append(current_hit)
                    else:
                        logger.info(f"    |--> The hit {str(current_hit)} has a coverage below the established threshold (coverage = {str(round(cov*100, 2))}).")
                # handle the case of no minimum coverage required (simply append hit)
                else:
                    if current_hit not in hits_dict[id]:
                        logger.info(f"There is no minimum coverage. Adding hit: {str(current_hit)}")
                        hits_dict[id].append(current_hit)
                
    logger.info(f"|--> Returning hits for {str(len(hits_dict))} records.")
    return hits_dict
