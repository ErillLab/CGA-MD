import subprocess
import tempfile
from Bio import SeqIO
import os

def launch_meme_search(configParams, output, meme_output = False):
    """
    Launch MEME motif discovery search.

    This procedure launches a MEME search using the specified parameters, and saves 
    the results to the specified output directory.

    Parameters:
        sequences_dict (dict): a dictionary containing the sequences to search for motifs.
            The keys represent sequence identifiers, and the values are BioPython SeqRecord objects.
        configParams (Variables): an instance of the Variables class containing MEME search parameters.
            The MEME search parameters include mod, nmotifs, minw, maxw, revcomp, and pal.
        output (str): the path to the output directory where MEME results will be saved.
        meme_output (bool or str): if False, output is the standard output.
            If not False, it indicates the output folder for the new MEME execution.

    Returns:
        list: a list of SeqRecord objects containing the input sequences in FASTA format.

    Note:
        The function generates a temporary FASTA file from the input sequences, launches the MEME search
        using the specified parameters, and saves the results to the specified output directory.
        If an error occurs during the process, appropriate error messages are printed.
    """
    
    #define location of output FASTA file
    output_file_path = os.path.join(output,"sequences_filtered.fasta")

    #set output folder for MEME if only rerunning MEME
    if meme_output: 
        output = os.path.join(output,meme_output)
    else:
        output = os.path.join(output,'meme_out')
    
    # Get MEME parameters
    mod = str(configParams.meme_mod)
    nmotifs = str(configParams.meme_nmotifs)
    minw = str(configParams.meme_minw)
    maxw = str(configParams.meme_maxw)
    options = ["-dna", "-mod", mod, "-nmotifs", nmotifs, "-minw", minw, "-maxw", maxw, "-brief", "10000"]
    
    #add additional parameters
    if configParams.meme_revcomp:
        options.append("-revcomp")
    
    if configParams.meme_pal:
        options.append("-pal")

    try:
        # Launch MEME search
        command = ["meme", output_file_path, "-oc", output] + options
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        
        # Check for errors
        if stderr:
            print("Error while running MEME:", stderr)
    except Exception as e:
        print("An error occurred:", e)
    