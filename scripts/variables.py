class Variables:
    """
    A class to represent variables used for storing config parameters.
    This class provides attributes to store various parameters and data used in the application.

    Attributes:
        REQUEST_LIMIT (int): the maximum number of requests allowed per unit of time for Entrez API.
        SLEEP_TIME (float): the time interval (in seconds) to sleep between consecutive requests to Entrez API.
        EMAIL (str): the email address associated with the Entrez API key.
        API_KEY (str): the API key required for accessing Entrez services.
        db (str): the name of the BLAST database to be used.
        tax_IDs (list): a list of taxonomic IDs for filtering BLAST search results.
        e_value (float): the E-value threshold for BLAST searches.
        query_coverage (float): the minimum coverage threshold for BLAST hits.
        max_hits (int): the maximum number of hits to be returned by a BLAST search.
        max_intergenic_size (int): the maximum size (in base pairs) allowed for intergenic regions.
        min_intergenic_size (int): the minimum size (in base pairs) allowed for intergenic regions.
        upstream_size_region (int): the size of the upstream region for sequence retrieval.
        downstream_size_region (int): the size of the downstream region for sequence retrieval.
        max_sequence_length (int): the maximum length of the sequence of the putative promoter.
        crawlback (bool): whether to crawl back to look for operon lead gene.
        max_identity (float): the maximum identity threshold for sequence filtering.
        meme_mod (str): the mod parameter for MEME search.
        meme_nmotifs (int): the number of motifs to search for in MEME.
        meme_minw (int): the minimum width of motifs to search for in MEME.
        meme_maxw (int): the maximum width of motifs to search for in MEME.
        meme_revcomp (str): whether to allow reverse complement motifs in MEME.
        meme_pal (str): whether to enforce sequence palindromes in MEME.
        run_only_meme (bool): indicates whether you want to run only meme or not
            (i.e., you already have the filtered sequences). If False, the program will be executed from the beginning.
        folder_name (str): the name of the output folder.
        input_records (list): a list containing input records or data.

    Note:
        Each attribute is initialized to None by default and should be assigned specific values as needed.
    """
    
    def __init__(self):
        # Entrez parameters
        self.REQUEST_LIMIT = None
        self.SLEEP_TIME = None
        self.EMAIL = None
        self.API_KEY = None
        # BLAST parameters
        self.db = None
        self.tax_IDs = None
        self.e_value = None
        self.query_coverage = None
        self.max_hits = None
        # Sequences parameters
        self.max_intergenic_size = None
        self.min_intergenic_size = None
        self.upstream_size_region = None
        self.downstream_size_region = None
        self.max_sequence_length = None
        self.crawl_back = None
        self.fuse_intergenic = None
        # Max identity
        self.max_identity = None
        # MEME parameters
        self.meme_mod = None
        self.meme_nmotifs = None
        self.meme_minw = None
        self.meme_maxw = None
        self.meme_revcomp = None
        self.meme_pal = None
        self.run_only_meme = None
        # Output parameters
        self.folder_name = None
        # Input records
        self.input_records = None
