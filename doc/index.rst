seq_crumbs
==========

Seq Crumbs aims to be a collection of small sequence processing utilities.

Seq Crumbs is modeled after the Unix command line text processing utilities so every utility tries to perform a specific task and most of them take a sequence file as input and create a new processed sequence file as output.
This design encorages the assembly of the Seq Crumbs utilities with Unix pipes.

-------

 
All Seq Crumbs try to share a consistent interface.
Most Seq Crumbs take can their input from the standard input to be able to work with Unix Pipes.
Alternatively several input sequence files can be provided as a list of arguments.
By default they throw their output to the standard output, although this behaviour can be changed with the *-o* parameter (or *--outfile*).

Seq Crumbs supports compressed gzip and BGZF_ files.
When used as input it autodetects the compressed files.
It can also generate compressed outputs.

The sequence formats accepted by Seq Crumbs are the ones supported by Biopython's SeqIO_ module.
As output only Sanger and Illumina fastq and fasta files are supported.

Seq Crumbs can take advantage of the multiprocessor computers by splitting the computational load into several processes.

------

 
Seq Crumbs is `free software`_. Licensed mainly under the General Public License (GPL_).
For more details on the licensing take a look at the LICENSE.txt file.

.. toctree::
   :maxdepth: 2
   :hidden:

   available_crumbs
   download
   install

.. include:: links.txt

