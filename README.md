# crAss-RPKM

A version of crAss (cross assembly) modified to support RPKM (Reads Per Kilobase per Million mapped reads).

This is particularly useful in viral metagenomics, and can be used to generate normalized heatmaps of contigs.


## Derived from crAss

This is derived from the original GPLv3 licensed code from:
Reference-independent comparative metagenomics using cross-assembly: crAss
Bas E. Dutilh et al. 2012

### Usage
Run crAss as follows:
perl crAss.pl input_directory

The input directory should contain the metagenomic reads (one file per metagenome)
and a cross-assembly ACE file. Reads files can be in FASTA format (extensions:
.fasta .fna or .fa) or FASTQ format (extensions: .fastq or .fq). Note that read
identifiers should be unique across all the files.

For further info please visit the crAss website at http://edwards.sdsu.edu/crass/


### Dependencies
The PERL script calls the following programs:
- BioNJ        http://www.atgc-montpellier.fr/bionj/
- GNUPlot      http://www.gnuplot.info/
- GhostScript  http://www.ghostscript.com/
- DrawTree     http://evolution.genetics.washington.edu/phylip/doc/drawtree.html
Make sure these programs are available and executable from your command line, in
the same directory as crAss.pl.


### License (GPLv3)
Copyright (C) 2012 Nick Semenkovich (@semenko) <semenko@alum.mit.edu>, Gordon Lab, Washington University School of Medicine in St. Louis
Copyright (C) 2012 Bas E. Dutilh

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
