from Bio import Entrez
import time
from logger import logger

def get_gene_record(hit, configParams):
    """"
    Retrieves the nucleotide record with the highest priority for a given list of protein accessions.

    Parameters:
        hit (str): accession number for the protein record of interest.
        configParams (Variables): an instance of the Variables class containing config parameters.

    Returns:
        max_p_record (dict or None): a dictionary containing information about the nucleotide record with the highest priority,
            or None if no gene records are found.

    Raises:
        IOError: if there's an issue with downloading the records.

    Notes:
        This function retrieves nucleotide records from NCBI based on protein accessions using Entrez.
        It then prioritizes the retrieved records based on predefined priority scores and returns
        the record with the highest priority.
    """
    scores = {
        "NC_": 7, "AC_": 7, "AE": 6, "CP": 6, "CY": 6,
        "NZ_": 5, "NT_": 5, "NW_": 5,
        "AAAA-AZZZ": 4, "AAAAAA-AZZZZZ": 4,
        "JAAA-JZZZ": 4, "JAAAAA-JZZZZZ": 4,
        "U": 3, "AF": 3, "AY": 3, "DQ": 3
    }

    genome_list = []

    logger.info(f"Retrieving nucleotide records for {hit}:")
    
    # obtain the IPG report for the protein ID
    for attempt in range(configParams.REQUEST_LIMIT + 1):
        try:
            handle = Entrez.efetch(db="protein", id=hit, rettype="ipg", retmode="xml")
            records = Entrez.read(handle)
            time.sleep(configParams.SLEEP_TIME)
            break
        except Exception as e:
            logger.info(f"|--> NCBI exception raised on attempt {attempt + 1}. Retrying.")
            if attempt == configParams.REQUEST_LIMIT:
                logger.info(f"Failed to download record after {configParams.REQUEST_LIMIT} attempts.")
                return None

    # assign priorities to genome records encoding the protein
    logger.info("|--> Recording gene records with priorities")
    if "IPGReport" in records and "ProteinList" in records["IPGReport"]:
        for protein in records["IPGReport"]["ProteinList"]:
            if "CDSList" in protein:
                for cds in protein["CDSList"]:
                    cds_acc = cds.attributes["accver"]
                    cds_start = cds.attributes["start"]
                    cds_stop = cds.attributes["stop"]
                    cds_strand = 1 if cds.attributes["strand"] == "+" else -1
                                            
                    cds_score = 0
                    for key, value in scores.items():
                        if cds_acc.startswith(key):
                            cds_score = value
                            break
                        
                    cds_record = {"acc": cds_acc, "start": cds_start, "stop": cds_stop,
                                  "strand": cds_strand, "score": cds_score}
                    genome_list.append(cds_record)
            else:
                continue
    else:
        logger.info(f"|--> No gene records found for {hit}")
        return None

    # return the "best" genome record encoding the protein
    if genome_list:
        max_p_record = max(genome_list, key=lambda x: x["score"])
        logger.info(f"|--> Returning gene {max_p_record["acc"]} record for {hit}. Priority score: {max_p_record['score']}")
        return max_p_record
    else:
        logger.info(f"|--> No gene records found for {hit}")
        return None
