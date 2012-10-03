=============
 sff_extract
=============

-----------------------------------
Extract reads from binary SFF files
-----------------------------------

:Author: jblanca@upv.es, pziarsolo@upv.es
:Date: 19-7-2012
:Version: 0.1
:Manual section: 1

SYNOPSIS
========

  **sff_extract** [**-o** *output_file*] [**--clip**] [**--min_left_clip**]
  [**--min_frequency**] *INPUT_SFFs* 


DESCRIPTION
===========

sff_extract extracts the reads from a list of SFF files and writes them in a sanger fastq file.


OPTIONS
=======

**-o**, **-output** *output_sequence_file*
		The path to the resulting output file, by default STDOUT.
		
**-c**, **--clip**
        Trim the reads by the SFF suggested quality and adaptor limits instead of masking them
        
**--min_left_clip** *number_of_base_pairs*
        Mask (or trim) at least the given number of nucleotides
        
**--max_percentage** *percentage*
        Maximum allowable percentage of reads starting with the same word before tagging the SFF file as "strange"
