re-pair
=======

Program to re-do the pairing of fastq reads. This program is modified from   http://code.google.com/p/ngopt/source/browse/trunk/tools/pair_reads/repair.cpp?r=85 


I added the ability to parse the new Illumina casava style headers. Also the program can handle R1/R3 as well as R1/R2.


Note that this program loads all reads into memory, so either run on a large memory system, or break down your fastq files into the smallest units possible where pairing is still guarenteed to exist.



