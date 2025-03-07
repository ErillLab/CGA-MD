from variables import Variables
from logger import logger

def extract_params(config):
    """
    Extract parameters from a config dictionary and store them in a Variables object.

    This function extracts specific parameters from the provided config dictionary and assigns them
    to corresponding attributes of a Variables object.

    Parameters:
        config (dict): a dictionary containing configuration parameters.

    Returns:
        Variables: an instance of the Variables class containing the extracted parameters.
    """
    
    configParams = Variables()
    
    # Entrez parameters
    configParams.REQUEST_LIMIT = config["entrez_parameters"]["request_limit"]
    configParams.SLEEP_TIME = config["entrez_parameters"]["sleep_time"]
    configParams.EMAIL = config["entrez_parameters"]["email"]
    configParams.API_KEY = config["entrez_parameters"]["api_key"]
    
    # BLAST parameters
    configParams.db = config["blast_parameters"]["db"]
    configParams.tax_IDs = config["blast_parameters"]["tax_IDs"]
    configParams.e_value = config["blast_parameters"]["e_value"]
    configParams.query_coverage = config["blast_parameters"]["query_coverage"]
    configParams.max_hits = config["blast_parameters"]["max_hits"]
    
    # Sequence parameters
    configParams.max_intergenic_size = config["sequences_parameters"]["max_intergenic_size"]
    configParams.min_intergenic_size = config["sequences_parameters"]["min_intergenic_size"]
    configParams.upstream_size_region = config["sequences_parameters"]["upstream_size_region"]
    configParams.downstream_size_region = config["sequences_parameters"]["downstream_size_region"]
    configParams.max_sequence_length = config["sequences_parameters"]["max_sequence_length"]
    configParams.crawl_back = config["sequences_parameters"]["crawl_back"]
    configParams.fuse_intergenic = config["sequences_parameters"]["fuse_intergenic"]
    #if not crawling back, no fusing possible
    if not configParams.crawl_back:
        print("|--> Crawl back disabled. Deactivating fuse intergenic option...")
        logger.info("|--> Crawl back disabled. Deactivating fuse intergenic option...")
        configParams.fuse_intergenic = False
    
    
    # Max identity (filter)
    configParams.max_identity = config["max_identity"]
    
    # MEME parameters
    configParams.meme_mod = config["meme_parameters"]["mod"]
    configParams.meme_nmotifs = config["meme_parameters"]["nmotifs"]
    configParams.meme_minw = config["meme_parameters"]["minw"]
    configParams.meme_maxw = config["meme_parameters"]["maxw"]
    configParams.meme_revcomp = config["meme_parameters"]["revcomp"]
    configParams.meme_pal = config["meme_parameters"]["pal"]
    configParams.run_only_meme = config["run_only_meme"]
    
    # Output parameters
    configParams.folder_name = config["output_parameters"]["folder_name"]
    
    # Input records
    configParams.input_records = config["input_records"]
    
    return configParams
