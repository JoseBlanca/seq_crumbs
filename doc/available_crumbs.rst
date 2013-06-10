Available Crumbs
----------------

sff_extract
    Extracts reads from an SFF_ file used by 454 and Ion Torrent.

split_matepairs
    Splits mate-pairs separated by an oligo sequence.

filter_by_quality
    Filters sequences according to mean quality.

filter_by_length
    Filters sequences according to maximum and minimum length thresholds.

filter_by_name
    Filters sequences with a list of names given in a file.

filter_by_blast
    Filters the sequences using BLAST.

filter_by_complexity
    Filters sequences according to their complexity.

trim_by_case
    Trims sequences according to case.

trim_edges
    Removes a fixed number of residues from sequence edges.

trim_quality
    Removes, using a sliding window, regions of low quality in the edges.

trim_blast_short
    Removes oligonucleotides by using the blast-short algorithm.

convert_format
    Converts between the different supported sequence formats.

guess_seq_format
    Guesses the format of a file, including Sanger and Illumina fastq formats.

cat_seqs
    Concatenates one or several input sequence files, possibly in different formats, into one output.

seq_head
    Outputs only the first sequences of the given input.

sample_seqs
    Outputs a random sampling of the input sequences.

change_case
    Modifies the case of sequences. Case can be converted to lower or upper, or swapped.

pair_matcher
    Filters out orphaned read pairs.

interleave_pairs
    Interleaves two ordered paired read files.

deinterleave_pairs
    Splits an ordered file of paired reads into two files, one for each end.

calculate_stats
    Generates basic statistics for the given sequence files.

orientate_transcripts
    Reverse complements transcripts according to polyA, ORF or BLAST hits.

fastqual_to_fastq
    Converts fasta and qual files to a fastq format file.
.. include:: links.txt
