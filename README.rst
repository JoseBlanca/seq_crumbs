
Seq Crumbs aims to be a collection of small sequence processing utilities.

Seq Crumbs is modeled after the Unix command line text processing utilities so every utility tries to perform a specific task and most of them take a sequence file as input and create a new processed sequence file as output.
This design encorages the assembly of the Seq Crumbs utilities with Unix pipes.

You can find more information about seq_crumbs in the seq_crumbs_ web site.


Available Crumbs
----------------

sff_extract
        It extracts the reads from a SFF_ file used for instance by 454 and Ion Torrent

split_matepairs
        It splits the mate-pairs separated by an oligo sequence.

filter_by_length
        It filters sequences according to a length threshold.

trim_by_case
        It trims the sequences according to its case.

trim_edges
        It removes a fixed number of residues from the sequence edges.

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

change_case
        It modifies the letter case of the sequences. The case can be converted to lower or upper case or swapped.

pair_matcher
        It filters out orphaned read pairs.

interleave_pairs
        It interleaves two ordered paired read files.

deinterleave_pairs
        It splits a ordered file of paired reads into two files, one for every end.


General Usage
---------------

All Seq Crumbs try to share a consistent interface.
Most Seq Crumbs take can their input from the standard input to be able to work with Unix Pipes.
Alternatively several input sequence files can be provided as a list of arguments.
By default they throw their output to the standard output, although this behaviour can be changed with the **-o** parameter (or **--outfile**).

Seq Crumbs compressed gzip and BGZF_ files.
When used as input they autodetect the compressed files and they can also generate compressed outputs.

The sequence formats accepted by Seq Crumbs are the ones supported by Biopython's SeqIO_ module.
As output only Sanger and Illumina fastq and fasta files are supported.

Seq Crumbs can take advantage of the multiprocessor computers by splitting the computational load into several processed.


Install
-------

Seq Crumbs depends on Python 2.7 and Biopython_.
The installation manual is located in the INSTALL document.


Related software
----------------

Seq Crumbs relies heavily on Biopython_ and without this free software project it won't be able to provide some of its functionalities.

Biopieces is project with a scope similar to Seq Crumbs.
Biopieces is a great software project capable of working with different kinds of biological data using Unix Pipes.
Seq Crumbs tries to be more limited in its scope limiting itself only to sequence files and thus providing a somewhat simpler interface.

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

