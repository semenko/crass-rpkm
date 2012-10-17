crAss
-----
Reference-independent comparative metagenomics using cross-assembly: crAss
Bas E. Dutilh et al. 2012


USAGE
-----
Run crAss as follows:
perl crAss.pl input_directory

The input directory should contain the metagenomic reads (one file per metagenome)
and a cross-assembly ACE file. Reads files can be in FASTA format (extensions:
.fasta .fna or .fa) or FASTQ format (extensions: .fastq or .fq). Note that read
identifiers should be unique across all the files.

For further info please visit the crAss website at http://edwards.sdsu.edu/crass/


DEPENDENCIES
------------
The PERL script calls the following programs:
- BioNJ        http://www.atgc-montpellier.fr/bionj/
- GNUPlot      http://www.gnuplot.info/
- GhostScript  http://www.ghostscript.com/
- DrawTree     http://evolution.genetics.washington.edu/phylip/doc/drawtree.html
Make sure these programs are available and executable from your command line, in
the same directory as crAss.pl.


COPYRIGHT AND LICENSE
---------------------
Copyright (C) 2012  Bas E. Dutilh

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
