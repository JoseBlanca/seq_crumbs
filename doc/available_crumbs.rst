Available Crumbs
----------------

sff_extract
  It extracts the reads from a SFF_ file used by 454 and Ion Torrent.

split_matepairs
    It splits the mate-pairs separated by an oligo sequence.

filter_by_quality
    It filters the sequences according to its mean quality.

filter_by_length
    It filters sequences according to maximum and minimum length thresholds.

filter_by_name
  It filters the sequences with a list of names given in a file.

filter_by_blast
    It filters the sequences using blast

filter_by_complexity
  It filters the sequences according to its complexity.

filter_by_bowtie2
  It filters the sequences using bowtie2

trim_by_case
    It trims the sequences according to its case.

trim_edges
    It removes a fixed number of residues from the sequence edges.

trim_quality
  It removes with a sliding windows regions of low quality in the edges

trim_blast_short
  It removes oligonucleotides by using the blast-short algoritm

convert_format
    It changes between the different supported sequence formats.

guess_seq_format
    It guesses the format of a file, including Sanger and Illumina fastq formats.

cat_seqs
    It concatenates one or several input sequence files, that can be in different formats, into one output.

seq_head
    It outputs only the first sequences of the given input.

sample_seqs
    It does a random sampling of the input sequences.

count_seqs
    It counts sequences in the input files

change_case
    It modifies the letter case of the sequences. The case can be converted to lower or upper case or swapped.

pair_matcher
    It filters out orphaned read pairs.

interleave_pairs
    It interleaves two ordered paired read files.

deinterleave_pairs
    It splits a ordered file of paired reads into two files, one for every end.

calculate_stats
    It generates basic statistics for the given sequence files

orientate_transcripts
    It reverse complements transcripts according to polyA, ORF or blast hits

fastqual_to_fastq
    It converts fasta and qual files to a fastq format file
.. include:: links.txt
