from logger import logger
from get_genome import fetch_genome_annotations
from Bio.SeqRecord import SeqRecord

def grab_intergenic_regions(query_id,hit, configParams):
    """
    Grabs intergenic regions in the downloaded genome fragment.

    This function identifies "valid" intergenic regions upstream of a gene hit based on the provided
    genome annotations and configuration parameters.

    Parameters:
        query_id: the id of the BLAST query that resulted in the hit
        hit (dict): a dictionary containing information about the gene hit, including start and stop positions,
            strand, and other relevant details.
        configParams (Variables): an instance of the Variables class containing configuration parameters,
            such as the minimum and maximum intergenic distances and region sizes.

    Returns:
        list: a list of annotated intergenic sequences (SeqRecord objects) upstream of the hit.
            
    Logic:
        The idea here is to identify the intergenic region(s) that is(are) most likely controling expression of
        the gene coding for an ortholog of the query protein.
        
        The most obvious candidate is the region directly upstream of the gene, but if the gene is preceded by
        other genes in the same direction (i.e. it is a directon), the idea is to crawl back until the likely
        operon promoter is identified.
        
        To cover all bases, the intergenic regions contained in the directon up to the "operon" lead gene are
        also grabbed.
        
        In this process, we use a min_intergenic_distance and a max_intergenic_distance.
        
        The min_intergenic_distance is the minimum intergenic distance that we are willing to consider as we
        crawl back. Any sequences under this value will be discarded (not grabbed).
        
        The max_intergenic_distance is the maximum intergenic distance that will determine that we have
        encountered the lead gene in an operon as we crawl back. If the size of the intergenic region between
        two consecutive genes in a directon is larger than this, the first gene is considered the lead operon gene.
    """
    
    # minimum intergenic distance in a directon required to grab intergenic sequence
    min_intergenic_distance = configParams.min_intergenic_size
    # maximum intergenic distance in a directon required to determine lead operon gene
    max_intergenic_distance = configParams.max_intergenic_size
    
    # maximum size for the lead intergenic region
    max_sequence_length = configParams.max_sequence_length
    
    # bool indicator for crawl back
    crawl_back = configParams.crawl_back

    # bool indicator for fusing together all intergenic regions for an "operon"
    fuse_intergenic = configParams.fuse_intergenic

    #list of intergenic regions grabbed
    grabbed_intergenic = []
    
    if hit["acc"] is None:
        logger.info("|--> No accession number provided in the gene hit.")
        return grabbed_intergenic
    
    logger.info(f"Grabbing intergenic regions from the genome of {hit['acc']}")
    
    #sequence to be analyzed
    seq = fetch_genome_annotations(hit, configParams)

    #get only gene features, and reverse the list so that ortholog-encoding gene is now first, not last
    gene_annotations = list(reversed([feature for feature in seq.features if feature.type == "gene"]))
    
    #go backwards through all gene features, grabbing intergenic regions as needed
    #until we hit operon lead gene (either d>max_intergenic or opposed strand)
    #or we go through the entire list of features
    for gene_index in range(len(gene_annotations)-1):
        #grab the region if it is larger than min intergenic distance
        #first get coordinates and compute size
        current_gene = gene_annotations[gene_index]
        next_gene = gene_annotations[gene_index + 1]
        #intergenic region size is distance from end of next to start of current
        int_start = next_gene.location.end
        int_end = current_gene.location.start
        int_region_size = int_end - int_start
        
        #grab protein ID, if available, for current gene, scanning CDS features
        current_gene.qualifiers['protein_id']=['None']
        for feature in seq.features:
            if feature.type == 'CDS':
                if 'locus_tag' in feature.qualifiers:
                    if feature.qualifiers['locus_tag'] == current_gene.qualifiers['locus_tag']:
                        if 'protein_id' in feature.qualifiers: #CDS records may not have protein_id (e.g. pseudogenes)
                            current_gene.qualifiers['protein_id'] = feature.qualifiers['protein_id']
                        #when we find a match betwen locus tags, we annotate the protein ID and break
                        break
        
        #if larger than min intergenic region
        if int_region_size > min_intergenic_distance:
            #check if region is larger than maximum allowed
            if int_region_size > max_sequence_length:
                #if so, set it to grab up to max_sequence_length upstream of current gene start
                int_start=int_end-max_sequence_length
            #create sequence record with intergenic region, labeled with protein and locus_tag
            curr_gene_locus_tag = current_gene.qualifiers['locus_tag'][0] if 'locus_tag' in current_gene.qualifiers else 'None'
            next_gene_locus_tag = next_gene.qualifiers['locus_tag'][0] if 'locus_tag' in next_gene.qualifiers else 'None'
            int_region = SeqRecord(seq.seq[int(int_start):int(int_end)],\
                                   id=query_id+'$'+current_gene.qualifiers['protein_id'][0]+'|'+ curr_gene_locus_tag + '@' + seq.id,\
                                   description=f"Region between {curr_gene_locus_tag} and {next_gene_locus_tag}",\
                                   name=current_gene.qualifiers['protein_id'][0])
                                   
            #append to list of intergenic regions
            grabbed_intergenic.append(int_region)
    
        #decide whether we stop loop (i.e. this was last intergenic region to grab)
        if (next_gene.strand == -1) or (int_region_size > max_intergenic_distance):
            if len(grabbed_intergenic)>0:
                logger.info(f"    |--> Intergenic regions found for {hit['acc']}: {len(grabbed_intergenic)}")
            break

        #break if we are only retrieving the first sequence (no crawl back)
        if not crawl_back:
            break

    #determine whether to return as independent sequences or as a "fused" sequence for the entire "operon"
    if fuse_intergenic:
        fi_seq = ''
        fi_id = query_id+'$'
        fi_desc = 'Region between ' + query_id+'$'

        #compose a combo sequence plus a summary header for the "fused" sequence        
        for inter_region in grabbed_intergenic:
            fi_seq = fi_seq + inter_region.seq
            fi_id = fi_id + inter_region.id.split('$')[1].split('|')[0] + '|' + str(len(inter_region.seq)) + '|'
            fi_desc = fi_desc + inter_region.id.split('$')[1].split('|')[0] + '|' + str(len(inter_region.seq)) + ', '
        fi_id = fi_id[:-1] + '@' + seq.id
        fi_desc = fi_desc[:-2] + ' in ' + seq.id
        
        grabbed_intergenic = [SeqRecord(fi_seq, id=fi_id, description=fi_desc)]
        
    return(grabbed_intergenic) 
    