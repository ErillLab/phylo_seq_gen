# phylo_seq_gen
Pulls sequences of a reference gene (i.e. recA) from a set of genomes. A MSA of these sequences can be used to construct a phylogeny.

## Usage

Need to setup the conda evironment.

Input file is expected to be in the ```./input/``` directory. Output is saved to the ```./output/``` directory. The output file is a FASTA file with the sequence of interest for each target species. 

```
python phylo_seq_gen.py [input_file_name] [output_file_name]
```

## Example Input file:

```
{
  "target_genomes": [
    {
      "Vibrio_cholerae_MS6_GCF_000829215.1": "GCF_000829215.1", 
      "Vibrio_vulnificus_YJ016_GCF_000009745.1": "GCF_000009745.1", 
      "Vibrio_metoecus_GCF_009665275.1": "GCF_009665275.1", 
      "Vibrio_vulnificus_GCF_002204915.1": "GCF_002204915.1", 
      "Vibrio_cidicii_GCF_001597945.1": "GCF_001597945.1"
    }
  ], 
  "query_file": "recA.fasta", 
  "api_key": "", 
  "email": ""
}
```
