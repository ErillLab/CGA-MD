from Bio import Entrez, SeqIO
import time
import logger

def fetch_genome_annotations(accession, configParams):
    """
    Retrieves the neighbouring genome sequence and annotations for a given accession number.

    Parameters:
        accession (dict): a dictionary containing the accession number under the key 'acc'.
        configParams (Variables): an instance of the Variables class containing configuration parameters.

    Returns:
        record (SeqRecord or None): a SeqRecord object containing the genome annotations if successful,
        or None if the retrieval fails.

    Raises:
        IOError: if there's an issue with downloading the annotations.

    Notes:
        This function retrieves the genome annotations from NCBI using Entrez based on the provided accession number.
        It will download the sequence record upstream_size_region bases upstream of the gene start and
        downstream_size_region bp downstream of the gene end.
        
        It donwloads the sequence based on the orientation of the gene we are interested in, so that the gene is ALWAYS
        the last gene of the returned sequence segment and it is ALWAYS in the forward strand (i.e. the NCBI efetch 
        command already returns the segment with the sequence (and annotations) already reverse complemented if needed.
        
        It attempts to download the annotations multiple times (up to REQUEST_LIMIT) if any errors occur.
        The retrieved annotations are returned as a SeqRecord object in GenBank format.
    
    """
    # read in the config parameters for strand, upstream and downstream distances to grab
    accession_strand = accession["strand"]
    upstream_size_region = configParams.upstream_size_region
    downstream_size_region = configParams.downstream_size_region


    # define the region to grab based on the strand we are targeting
    # adding the appropriate upstream/downstream extensions
    if accession_strand == 1:
        strand = 1
        start = int(accession["start"])-upstream_size_region
        end = int(accession["stop"])+downstream_size_region
    elif accession_strand == -1:
        strand = 2  #<-- this will instruct efetch to return the reverse-complement of the sequence segment
        start = int(accession["start"])-downstream_size_region
        end = int(accession["stop"])+upstream_size_region

        
    for attempt in range(configParams.REQUEST_LIMIT + 1):
        try:
            handle = Entrez.efetch(db="nuccore", id=accession["acc"], seq_start=start, seq_stop=end, strand=strand, rettype="gbwithparts", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            time.sleep(configParams.SLEEP_TIME)
            break
        except Exception as e:
            logger.info(f"NCBI exception raised on attempt {attempt + 1}. Retrying.")
            if attempt == configParams.REQUEST_LIMIT:
                logger.info(f"Failed to download record after {configParams.REQUEST_LIMIT} attempts.")
                return None
    
    return record
