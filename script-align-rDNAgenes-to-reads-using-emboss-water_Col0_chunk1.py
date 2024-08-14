import gzip
from Bio import SeqIO
import os
import pandas as pd
import subprocess
from script_handle_files import (
    extract_file_name_without_suffix,
    extract_sequence_names,
    init_csv,
    append_csv,
    set_nan_values,
)

def summary_genes(file, output):
    """
    Summarize genes from a FASTA file and write to a CSV file.
    Args:
        file (str): Path to the input FASTA file.
        output (str): Path to the output CSV file.
    """
    summary_genes = []
    with open(file, "rt") as handle:
        for i,record in enumerate(SeqIO.parse(handle, "fasta")):
            summary_genes.append([record.name, len(record.seq)])

    results_data_frame = pd.DataFrame(summary_genes)
    # Column names
    columns = ["gene", "seq_length"]
    # Setting column names after DataFrame creation
    results_data_frame.columns = columns

    # Writing the .csv file
    try:
        results_data_frame.to_csv(output)
        print(f"CSV file '{output}' written successfully.")
    except Exception as e:
        # If an exception occurs during the try block, catch it and raise a new exception
        raise Exception(f"Error writing CSV file '{output}': {e}")

def run_command(command):
    """
    Run a shell command.
    Args:
        command (str): Command to execute.
    """
    try:
        # Use subprocess.run to execute the command
        result = subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        # Handle exceptions for non-zero exit codes
        print(f"Error running command: {e}")

def run_water(aseq, bseq, output):
    """
    Run water method from emboss software.
    Args:
        aseq (str): Path to the A sequence file.
        bseq (str): Path to the B sequence file.
        output (str): Path to the output file.
    """
    method = "water"
    #set considered sequence files
    arguments = []
    option_input = f"-asequence {str(aseq)} -bsequence  {str(bseq)}"
    arguments.append(option_input)
    
    #set methods parameter
    gapopen_value = 25
    gapextend_value = 1
    option_values=f"-gapopen {str(gapopen_value)} -gapextend {str(gapextend_value)}"
    arguments.append(option_values)
    
    #set output file
    option_output = f"-outfile {output}"
    arguments.append(option_output)

    #construct command
    water_command =f"{method} {' '.join(arguments)}"
    #print(water_command)
    run_command(water_command)
    
def extract_results_for_gene(file_path, gene_name):
    """
    Extract alignment results for a specific gene.
    Args:
        file_path (str): Path to the alignment file.
        gene_name (str): Name of the gene.
    Returns:
        list: List containing alignment information for the gene.
    """

    alignment_info = []
    line_start = f"# 2: {gene_name}" 
    start_idx = -1
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if line.startswith(line_start):
            alignment_info.append(gene_name)
            start_idx = i + 1
            break
    
    if start_idx == -1:
            raise ValueError("Gene name not found in file.")
            
    # Get the lines in between the start and end markers
    info = lines[start_idx : start_idx + 11]

    for line in info:
        if line.startswith("# Length:"):
            alignment_info.append(line.split(":")[1].strip())
        elif line.startswith("# Identity:"): 
            identity_info = line.split(":")[1].strip().split("/")
            alignment_info.append(identity_info[0].strip())
        elif line.startswith("# Similarity:"):
            similarity_info = line.split(":")[1].strip().split("/")
            alignment_info.append(similarity_info[0].strip())
        elif line.startswith("# Gaps:"):
            gaps_info = line.split(":")[1].strip().split("/")
            alignment_info.append(gaps_info[0].strip())
        elif line.startswith("# Score:"):
            alignment_info.append(line.split(":")[1].strip())

    return alignment_info



def extract_results_alignment(alignment_file, gene_names):
    """
    Extract alignment results for multiple genes.
    Args:
        alignment_file (str): Path to the alignment file.
        gene_names (list): List of gene names.
    Returns:
        list: List containing alignment information for all genes.
    """
    result = []
    for gene in gene_names:
        alignment_info = extract_results_for_gene(alignment_file, gene)
        result.extend(alignment_info)
    return result



def  main(fastq_gz_file_path, file_genes, output_dir):
    """Main program"""
    # Define the number of sequences to read
    num_sequences_to_read = 10
    
    # Define minimum sequence length 
    min_seq_len = 5000

    # extract data name and determine the output files    
    data_name = extract_file_name_without_suffix(fastq_gz_file_path)
    help_output = f"{output_dir}{data_name}_read_sequence.fastq"

    # extract gene names
    gene_names = extract_sequence_names(file_genes)

    #initialize csv file
    result_csvfile=f"{output_dir}{data_name}_alignment_info.csv"
    init_csv(result_csvfile,gene_names)

    # Open the compressed FASTQ file and parse its records using BioPython
    with gzip.open(fastq_gz_file_path, "rt") as handle:
        for i,record in enumerate(SeqIO.parse(handle, "fastq")):
            if i > num_sequences_to_read :
                # Open help output file and write read sequence into it
                with open(help_output, "w") as output_handle:
                    SeqIO.write(record, output_handle, "fastq")

                #initialize list for results
                results = []
                # Access information about each sequence
                results.append(str(record.name))
                seq_len = len(record.seq)
                results.append(seq_len)

                if seq_len >= min_seq_len: 
                    output_file=f"{output_dir}{data_name}_alignment_results.txt"
                    run_water(help_output, file_genes, output_file)
                    #time.sleep(20)
                    alignment_result = extract_results_alignment(output_file, gene_names)
                    #print(alignment_result)
                else:
                    alignment_result = set_nan_values(gene_names)

                results.extend(alignment_result)
                append_csv(result_csvfile, results)

            
            # # Break out of the loop when num_sequences_to_read is reached
            # if i + 1 == num_sequences_to_read:
            #     break





if __name__ == "__main__":
    
    """ Set directories and input files """
    # Get the current working directory
    current_directory = os.getcwd()

    # Print or use the current working directory
    #print("Current Working Directory:", current_directory)

    # Specify the path to your compressed FASTQ file including the considered reads
    input_dir = os.path.join(current_directory, "../../", "data/")
    input_fastq_gz_file_path = os.path.join(input_dir, "Col0_split_chunk_1.fastq.gz")
  
    # Specify the genes
    genes_dir = os.path.join(input_dir, "rDNA_unit_genes/")
    file_genes=f"{genes_dir}/all_genes_seqs_renamed.fasta"
    #print(file_genes)

    # Specify the path to output FASTQ file
    output_dir = os.path.join(current_directory, "../../results/identify_rDNA/")
    os.makedirs(output_dir, exist_ok=True)
    summary_genes_file= f"{output_dir}Col0_chunk1.csv"

    summary_genes(file_genes, summary_genes_file)
    main(input_fastq_gz_file_path, file_genes, output_dir)


