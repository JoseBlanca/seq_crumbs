
Seq Crumbs is a collection of small sequence processing utilities.

Seq Crumbs is modeled after the Unix command line text processing
utilities, so every utility tries to perform a specific task. Most of them
take a sequence file as input and create a new processed sequence file as
output.  This design encourages chaining the operation of multiple Seq
Crumbs utilities using Unix pipes.

You can find more information about seq_crumbs in the seq_crumbs_ web site.


Available Crumbs
----------------

sff_extract
    Extracts the reads from a SFF_ file used by 454 and Ion Torrent.

split_matepairs
    Splits the mate-pairs separated by an oligo sequence.

filter_by_quality
    Filters the sequences according to its mean quality.

filter_by_length
    Filters sequences according to a length threshold.

filter_by_name
    Filters the sequences with a list of names given in a file.

filter_by_blast
    Filters the sequences using BLAST.

filter_by_bowtie2
    Filters the sequences using bowtie2.

filter_by_complexity
    Filters the sequences according to their complexity.

trim_by_case
    Trims the sequences according to case.

trim_edges
    Removes a fixed number of residues from the sequence edges.

trim_quality
    Removes, using a sliding window, regions of low quality in the edges.

trim_blast_short
    Removes oligonucleotides by using the blast-short algorithm.

convert_format
    Converts between different supported sequence formats.

guess_seq_format
    Guesses the format of a file, including Sanger and Illumina fastq formats.

cat_seqs
    Concatenates one or more input sequence files, possibly in different formats, into one output.

seq_head
    Outputs only the first sequences of the given input.

sample_seqs
    Outputs a random sampling of the input sequences.

change_case
    Modifies the letter case of the sequences. The case can be converted to lower or upper case or swapped.

pair_matcher
    Filters out orphaned read pairs.

interleave_pairs
    Interleaves two ordered paired read files.

deinterleave_pairs
    Splits an ordered file of paired reads into two files, one for each end.

calculate_stats
    Generates basic statistics for the given sequence files.

count_seqs
    Counts the number of sequences and the total sequence length for the given files.

orientate_transcripts
    Reverse complements transcripts according to polyA, ORF or blast hits.

General Usage
---------------

All seq crumbs try to share a consistent interface.  By default most Seq
Crumbs read from standard input and write to standard output, allowing them
to to be easily combined using Unix pipes.  Alternatively, several input
sequence files can be provided as a list of arguments.  Output can also be
directed to specific files with the *-o* parameter (or *--outfile*).

seq_crumbs supports compressed gzip, BGZF_ and bzip2 files.
When used as input it autodetects the compressed files.
It can also generate compressed outputs.

The sequence formats accepted by seq_crumbs are those supported by Biopython's SeqIO_ module.
As output only Sanger and Illumina fastq and fasta files are supported.

seq_crumbs can take advantage of multiprocessor computers by splitting the computational load into several processes.

The filtering seq crumbs can be made aware of paired reads and can filter both reads of pairs at once.

Installation
------------

seq_crumbs depends on Python 2.7. Biopython_ is a recommended dependency.
The installation manual is located in the doc/install.rst document.


Related software
----------------

seq_crumbs relies heavily on Biopython_ and without this free software project it wouldn't be able to provide some of its functionality.

Biopieces_ is a project with a scope similar to seq_crumbs.
Biopieces_ is a great software project for working with different kinds of biological data using Unix pipes.
seq_crumbs tries to be more limited in its scope, limiting itself only to sequence files and thus providing a somewhat simpler interface.

Another software package very similar in approach to seq_crumbs is the nice fastx_ collection.

Other related software: PRINSEQ_, ea-utils_, Pyrocleaner_, `Sequence Cleaner <http://seqclean.sourceforge.net/>`_, lucy_, `NGS QC Toolkit <http://www.nipgr.res.in/ngsqctoolkit.html>`_, scythe_, sickle_, cutadapt_, trimomatic_ and FastQC_.

License
-------

Seq Crumbs is `free software`_. Licensed mainly under the General Public License (GPL_).
For more details on the licensing take a look at the LICENSE.txt file.


.. _seq_crumbs: http://bioinf.comav.upv.es/seq_crumbs/
.. _SFF: http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=format#sff
.. _BGZF: http://samtools.sourceforge.net/SAM1.pdf
.. _SeqIO: http://biopython.org/wiki/SeqIO
.. _Biopython: http://biopython.org/wiki/Biopython
.. _free software: http://en.wikipedia.org/wiki/Free_software
.. _GPL: http://www.gnu.org/copyleft/gpl.html
.. _fastx: http://hannonlab.cshl.edu/fastx_toolkit/
.. _Biopieces: http://code.google.com/p/biopieces/
.. _Pyrocleaner: https://pyrocleaner.mulcyber.toulouse.inra.fr/plugins/mediawiki/wiki/pyrocleaner/index.php/Pyrocleaner
.. _lucy: http://lucy.sourceforge.net/
.. _FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _scythe: https://github.com/vsbuffalo/scythe
.. _sickle: https://github.com/najoshi/sickle
.. _cutadapt: http://code.google.com/p/cutadapt/
.. _PRINSEQ: http://prinseq.sourceforge.net/manual.html
.. _trimomatic: http://www.usadellab.org/cms/index.php?page=trimmomatic
.. _ea-utils: http://code.google.com/p/ea-utils/
