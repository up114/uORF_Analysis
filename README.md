# uORF_Analysis

**all_uorfs_FINAL.py**
Data files required
- Transcript file:
    - in fa.gz file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz*

**count_reads.py**
Data files:
- A .ribo ribosome profiling file with coverage data
    - Example .ribo files found in **Data Files**
      
**filter_uorfs_FINAL.R**

**check_periodicity.py**
**format_uorfs.py**
Data files:
- Transcript file:
    - in fa.gz file format
    - Genetic sequences should be preceded with gene information formatted as: *>**ENSMUST00000070533.4**|ENSMUSG00000051951.5|OTTMUSG00000026353.2|OTTMUST00000065166.1|**Xkr4-201**|Xkr4|3634|**UTR5:1-150**|**CDS:151-2094**|UTR3:2095-3634|*
    - Sequences in CAPITAL LETTERS
    - like *appris_mouse_v2_selected.fa.gz*
      
**transcript2genome.py**
Data files:
- GENCODE GTF file:
    - Transcript names must match the version in the transcript file
    - for *appris_mouse_v2_selected.fa.gz,* use GENCODE vM25 (found at https://www.gencodegenes.org/mouse/release_M25.html)

**uorf_reads_psite.py**
Data files"
- Studies csv file for RiboBase:
    - should include study name, experiment name, and cell line columns
    - like *mouse_filtered_complete.csv*

