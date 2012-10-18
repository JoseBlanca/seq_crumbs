
seq_crumbs
==========

seq_crumbs aims to be a collection of small sequence processing utilities.

seq_crumbs is modeled after the Unix command line text processing utilities so every utility tries to perform a specific task and most of them take a sequence file as input and create a new processed sequence file as output.
This design encourages the assembly of the seq_crumbs utilities with Unix pipes to create complex pipelines.
There are already other similar software that uses this approach like `fastx <http://hannonlab.cshl.edu/fastx_toolkit/>`_ and `biopieces <http://code.google.com/p/biopieces/>`_.

seq_crumbs is a young software, that we use in our analysis and that we share with the hope that could be useful to you.
We have tried to create automated tests for every feature, but we are sure that some hidden bugs remain, so if you use it and you come across some problem we would appreciate your report.
You can also track or participate in the development of seq_crumbs in `github <https://github.com/JoseBlanca/seq_crumbs>`_.

 
All seq_crumbs try to share a consistent interface and implementation:
  * All take a sequence file as the input and most generate a processed sequence file as output.
  * The input can be either a file or the standard input and the output can also be a file or the standard output.
  * They can autodetect and work with files compressed with either gzip or BGZF_
  * The can work with any format supported by biopython's SeqIO_ and they try to autodetect the most common formats: fasta, Sanger and Illumina fastq.
  * Most seq_crumbs can split the work load in multicore machines into several processes

 
seq_crumbs is powered by Biopython_ library.

seq_crumbs is `free software`_. Licensed mainly under the General Public License (GPL_).
For more details on the licensing take a look at the LICENSE.txt file included in the seq_crumbs distribution.

.. toctree::
   :maxdepth: 2
   :hidden:

   available_crumbs
   download
   install

.. include:: links.txt

