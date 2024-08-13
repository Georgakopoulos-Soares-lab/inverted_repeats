# Inverted Repeats Analysis


# Extraction

For the extraction we used non-B gfa tool.

```
https://github.com/abcsFrederick/non-B_gfa
```

In particular, a python wrapper was built around non-B gfa, which post-processed the output file into a more informed tabular format.

# How to use the wrapper file

After cloning the repository, inside the scripts/ there is the wrapper file. To quickly run the snakemake workflow, you can follow the following procedure.  
cd into the scripts/ directory and run configure.sh after you grant the appropriate permissions (e.g. chmod 777). This will install the dependencies for the project.
Bear in mind I have not tested this on a different computer, but only for the purposes of the research. 

After all the required packages have been installed, you need a ```filtered_accessions.txt``` file, which has all the FULL paths of genomic assemblies in fasta format, on which we will extract the IR sequences. Also, inside the .env file you need to pass to the variable nonBDNA to the path of the gfa executable that we cloned from the previous step. Now, we can run the pipeline, as follows:

```
snakemake --snakefile nonbdna_pipe.py --cores 1
```

This will output the IR sequences, in the following directory:

```
inverted_repeats/scripts/repeat_out
├── IR_completed
│   ├── bucket_0.IR.completed
│   ├── bucket_1.IR.completed
│   ├── empty_accessions.IR.txt
│   └── inverted_repeats_db.parquet.snappy
├── IR_extracted_accessions
│   └── GCF_000002515.2_ASM251v1_genomic.IR.csv
├── IR_temp
│   └── GCF_000002515.2_ASM251v1_genomic_IR.tsv
└── schedule_2.json

```


# Analysis

The analysis was conducted on four stages:

- Basic EDA analysis, including densities across domains, kingdoms and phylums.
- Most prevalent inverted repeat arms amongst the four domains, using wordclouds, venn diagrams and upset plots.
- Biophysical Properties of IR sequences. The GC content and dinucleotide composition of all IR arms was examined and compared to the base ratio. Another analysis included the arm length and spacer length of IRs with varying arm length. 
- Compartment analysis. A snakemake workflow pipeline was built in order to extract the coverage of IRs with genomic subcompartments of interest, such as 
CDS, exons, genes.
- TSS/TES positioning of IRs. A pipeline was engineered in order to extract the relative positioning of IRs in a 1kb window around TSS/TES. The coordinates were extracted and processed from the publicly available GFF files.

# Contact

For any questions or bugs, please contact:

```
Nikol Chantzi; nmc6088@psu.edu
Ilias Georgakopoulos-Soares; izg5139@psu.edu
```
