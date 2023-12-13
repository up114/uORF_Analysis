# uORF_Analysis
 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**Environment**
- ribo_env.yml: contains the necessary packages for all function other than the genomic overlap comparison
- bedtools_env.yml: for the genomic overlap comparison (contains bedtools)
 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**all_uorfs_FINAL.py**

Data files required
- Transcript file:
    - in Fasta GZ file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz*

Notes:
The file with all of the uORF candidates can be found in Data Files > ltdstart_uorf.csv.

ltdstart_uorfs.csv: 
-   Mouse APPRIS transcript
-	Start: ["ATG", "ACG", "ATT", "CTG", "GTG", "TTG"]
-	Stop: ["TAA", "TGA", "TAG"]
-	Length range: 33-300 nucleotides (10-99 AA)
  
If you want to change the start codons, stop codons, or length range, you will need to run Code > all_uorfs_FINAL.py. Edit the start and stop codons in find_uORFS(), and edit the length when write_uORFS() is called at the bottom of the program. Make sure to have the correct path to the APPRIS Fasta GZ file – the APPRIS file can be found in Data Files. 

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**count_reads.py**

Data files:
- A .ribo ribosome profiling file with coverage data
    - Example .ribo files found in **Data Files**

Notes: 
This is the program for compiling the reads for all of the uORFs, and the code to run it is found at the bottom of the program file. Make sure to update the ribo_path variable. To run the function, you have two arguments: 
1)	File – this is ltdstart_uorf_no_overlap.csv or the outputted file from all_uorfs_FINAL.py
2)	Outfile – this is the name of the outputted csv file – e.g. “neural_reads.csv”

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**filter_uorfs_FINAL.R**

Notes: This program filters the reads results by CPM and start codon (with an ATG preference), and this program is run in _RStudio_. The You will have two arguments:

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**check_periodicity.py**

Notes: This program determines the periodicity output, and the code to run it is found at the bottom of the program file.You have three arguments:
1) The csv file from the preceding R Script 
2) The name of the complete csv output file (periodicity results for all inputted uORFs)
3) The accepted csv output file (periodicity results for the uORFs with p < 0.05).

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**format_uorfs.py**

Data files:
- Transcript file:
    - in Fasta GZ file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz*

Notes: This program adds additional information to an file containing uORFs. The added informations is: the transcript name, the Kozak sequence, the nucleotide sequence, the complete sequence (the sequence from the start of the transcript until right before the stop codon), the amino acid sequence, and the CDS start and stop sites.

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**transcript2genome.py**

Data files:
- GENCODE GTF file:
    - Transcript names must match the version in the transcript file
    - for *appris_mouse_v2_selected.fa.gz,* use GENCODE vM25 (found at https://www.gencodegenes.org/mouse/release_M25.html)
- Transcript file:
    - in Fasta GZ file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz*
- input sequences file:
    - in csv format
    - The input transcript pos table is seperated by comma. The input sequences file does not need to be edited if run through the rest of the pipeline. An example format of input file:
       - Gene,Start,Transcript,Stop
       - 1110012L19Rik,434,ENSMUST00000053981.5,434
       - 1110012L19Rik,768,ENSMUST00000053981.5,768
       - 1110012L19Rik,981,ENSMUST00000053981.5,981
         
Notes: This script converts the uORF transcript coordinates to genomic coordinates. 

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**Running the genomic overlap comparison**

Data files:
- GENCODE GTF file:
    - Transcript names must match the version in the transcript file
    - for *appris_mouse_v2_selected.fa.gz,* use GENCODE vM25 (found at https://www.gencodegenes.org/mouse/release_M25.html)
    - **Please note that you need to run the original GTF file through the make_cds_file() function in transcript2genome.py to filter the GTF file to only contain CDS entries**
- Transcript file:
    - in Fasta GZ file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz*
      
Notes: The genomic comparison only needs to be run on the nonoverlapping uORFs, so filter the uORFs before running the following steps. After converting the uORFs to genomic coordinates, you will need to convert the uORF genomic coordinates into a BED formatted file. I used Notepad++ to do so, and an example BED file can be found in the Data Files > ltdstart_bed.bed. 
Finally, overlap between GENCODE CDS regions and the identified uORFs can be found using bedtools. Run the make_cds_file() function in transcript2genome.py to filter the GTF file to only contain CDS entries. The instructions are as follows:

- Activate a Conda environment with bedtools - see bedtools_env.yml: 
     - conda activate bedtools_env
Run the following in Linux (I used Ubuntu):
- To get rows of overlap from the Gencode GTF: 
   -	 bedtools intersect -u -a "_path_/gencode_cds.gtf" -b "_path_/ltdstart_bed.bed" > "_path_/gtf_overlap.txt"
- To get rows of overlap from the uORF BED: 
   -	 bedtools intersect -u -a "_path_/ltdstart_bed.bed" -b "_path_/gencode_cds.gtf" > "_path_/uorf_overlap.txt"    

Note: These functions will report if at least one overlap is found. To also report the number of bp that overlap, add "-wo" to the options

 ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
**uorf_reads_psite.py**

Data files:
- Studies CSV file for RiboBase
    - should include study name, experiment name, and cell line columns - make sure to edit the specific column names in the script
    - like *mouse_filtered_complete.csv*
 - uORF CSV file
    - you will need to reformat your file to only contain genes, and start & stop transcript coordinates (to avoid an extremely large output). Example formatted file:
       - Gene,Start,Stop
       - 1110012L19Rik,434,512
       - 1110012L19Rik,101,113
       - 1110012L19Rik,981,1002

Notes: You will need to edit the file to include the correct paths in your TACC environment.  To edit, check all places where file paths are mentioned in the uorf_reads() function; I have put #EDIT comments next to these lines. The function called is located at the bottom of the scripts and requires four inputs:
1)	Database: The csv containing the studies to be examined
2)	Transcripts: a csv file with the transcript coordinate of the uORFs; find an example of the format of the file in the Data Files > ltdstart_gene_codons.csv (this csv should have no header)
3)	Outfile: the path for the output csv file with the Ribobase reads
4)	Outfile_offset: the path for the output csv file with the p-site offsets for each experiment

To run: Edit the sbatch script in Code > sbatch.sh to include the correct path for your uorf_reads_psite.py file and the correct conda environment. Then, run the job in TACC (I use LS6), using the command: sbatch sbatch.sh. 
