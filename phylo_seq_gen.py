'''
@author: ichaudr

Objective: Pull sequences of a reference gene from a set of input organisms that can be used to generate a phylogeny.
Note: This script is standalone from the operon_conserve_detect script; however, the output from that script can be reformatted to be the input 
    for this script. 

Input file: JSON

Workflow:
    Import the target genomes information from the JSON file.
    Download the Genebank records for those genomes and store them locally. 
    Assemle a BLAST db using the local Genebank records.
    Conduct BLAST searchs against the local DB restricted against the TAXID of each of the target genomes. 
    Retrieve the nucelotide sequence that got the hit, and store to output CSV file. 

'''

from Bio import Entrez, SeqIO, SeqRecord
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from tqdm import tqdm
from datetime import datetime
import os
import json
import time
import csv
import sys

#Input JSON file
input_file = ''

#Output CSV file
output_file = ''
output_records = []

#Dict of the target genomes where {key:value} -> {"species_name":"genome_assembly_accession"}
target_genomes = {}

#Local directory and file name template to download the Genbank files for the target genomes. 
cache_file = './cache/{species}.fasta'

#Local directory and template for BLAST databases
db_file = './b_db/{species}/{species}_blast-db'

#Entrez parameters
REQUEST_LIMIT = 5
SLEEP_TIME = .5
Entrez.email = ''
Entrez.api_key = ''

#BLAST variables
query_protein_path = './query/recA.fasta'


def blast_search(sp):
    '''
    Conducts the BLAST search for the specified species

    Parameters
    ----------
    sp: string
        Name of the species

    Returns
    -------
    None
    '''
    
    #Template command
    cmd = 'blastp -query {q_path} -out {out_path} -db {db_path} -outfmt 5'

    #Set up the parameters
    temp_out_path = './temp'
    db_path = '"' + db_file.format(species=sp) + '"'

    #Conduct the search
    tqdm.write('\tBLAST: ' + str(sp) + '...')
    os.system(cmd.format(q_path=query_protein_path, out_path=temp_out_path, db_path=db_path))

    #Storing the final results for this species
    tqdm.write('\tProcessing BLAST results ' + str(sp) + '...')

    with open(temp_out_path) as file:

        #Parse the output file from the BLAST search
        result_handle = list(NCBIXML.parse(file))
        
        if len(result_handle[0].alignments) > 0:
            record = result_handle[0].alignments[0]
        
            #Extract the accession for the alignment
            target_id = record.hit_def

            #Get the exact record that corresponds to this alignment accession
            seq_records = list(SeqIO.parse(cache_file.format(species=sp), 'fasta'))

            #Keeps track of whether or not a sequence has been added for this species
            added = False

            for seq_record in seq_records:
                if seq_record.description == target_id and added == False:
                    #output_description = '|' + '_'.join(str(e) for e in ((((target_id.split(' ')[0]).split('|')[1]).split('_'))[0:2]))
                    tqdm.write('\tAppending output record with ' + str(sp))
                    output_records.append(SeqRecord.SeqRecord(seq_record.seq, id=sp, description='|'))
                    added = True


def compile_database(sp):
    '''
    Assembles the downloaded target genomes for the requested species into one BLAST DB per species. 

    Parameters
    ----------
    sp: string
        Name of the species

    Returns
    -------
    None
    '''

    #Command template
    cmd = 'makeblastdb -in {input_file} -out {output_file} -dbtype prot'
    
    tqdm.write('\tCompiling database for' + str(sp) + '...')

    #Set up formatted file directories
    input_file = '"' + cache_file.format(species=sp) + '"'
    output_file = '"' + db_file.format(species=sp) + '"'

    #Execute command
    os.system(cmd.format(input_file=input_file, output_file=output_file))


def download_target_genomes(sp):
    '''
    Downloads the target genomes so they can be added to a local database for a local BLAST search.

    Parameters
    ----------
    sp: string
        Name of the species

    Returns
    -------
    None
    '''

    #Checks to see if this record has already been downloaded. 
    if os.path.exists(cache_file.format(species=sp)):
        tqdm.write(str(sp) + " already downloaded, continuing...")
        return

    #Sets the out
    output_file = cache_file.format(species=sp)

    #Gets all genomes associated with the species name by searching the nucleotide DB with the genome assembly accession
    for i in range(REQUEST_LIMIT):

        try:
            handle = Entrez.esearch(db='nuccore', term=(target_genomes[sp]), retmode='xml', retmax='5000')
            search_records = Entrez.read(handle)
            time.sleep(SLEEP_TIME)
        except:
            print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

            if i == ( REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str( REQUEST_LIMIT) + " attempts")
    


    #Fetch the accession from the search ID
    tqdm.write("\tDownloading " + str(sp) + "...")

    #Stores all the sequence objected for all the proteins in the species
    sequences_to_write = []

    #Fetches the records
    for i in range(REQUEST_LIMIT):

        try:
            handle = Entrez.efetch(db='nuccore', id=search_records['IdList'], start='begin', stop='end', rettype='fasta_cds_aa')
            records = list(SeqIO.parse(handle, 'fasta'))
            time.sleep(SLEEP_TIME)
            break
        except:
            print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now for ...")

            if i == ( REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str( REQUEST_LIMIT) + " attempts")

    tqdm.write('\tExtracting records for ' + str(sp))
    sequences_to_write.extend(records)
    
    if len(sequences_to_write) > 0:
        #Write all the sequences to the output file
        tqdm.write('\tWriting all CDS features to file for ' + str(sp))
        SeqIO.write(sequences_to_write, output_file, 'fasta')
    else:
        tqdm.write('\tNo sequences found for ' + str(sp))


def load_input_file():
    '''
    Loads all information from the input file assigned above.

    Parameters
    ----------
    None

    Returns
    -------
    None
    '''

    #File reader to extract data from the JSON
    file_reader = json.load(open('./input/' + input_file))

    #Load the target genome
    global target_genomes
    target_genomes = file_reader['target_genomes'][0]

    #Load the query file
    global query_protein_path
    query_protein_path = './query/' + file_reader['query_file']

    #Load Entrez parameters
    Entrez.api_key = file_reader['api_key']
    Entrez.email = file_reader['email']
    


if __name__ == "__main__":

    #Take note of the start time
    start_time = datetime.now()

    #Parse out the input filename and output file name
    if len(sys.argv) >= 3:        
        input_file = sys.argv[1]        
        output_file = sys.argv[2]
    else:
        input_file = 'rvva_input.json'
        output_file = 'output_rvva_10072020_2.fasta'


    #Load input file
    load_input_file()

    #Iterate through all the species and process them
    for sp in tqdm(target_genomes.keys(), desc='Processing'):
        tqdm.write(str(sp) + '\n' + '-'*50)
        download_target_genomes(sp=sp)

        if os.path.exists(cache_file.format(species=sp)):
            compile_database(sp=sp)
            blast_search(sp=sp)

        tqdm.write('\n')

    #Write to file
    SeqIO.write(output_records, ('./output/' + output_file), 'fasta')

    #Take note of the end time
    end_time = datetime.now()

    #Calculate elapsed time and print out
    elapsed_time = end_time - start_time
    print("*"*10 + " Finished in " + str(elapsed_time) + "*"*10)
