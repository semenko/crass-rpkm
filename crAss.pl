#!/usr/bin/perl -w
use strict;
use Math::Trig;
use Data::Dumper;

# crAss version 1.1-RPKM
# (c) 2012, Nicholas P. Semenkovich <semenko@alum.mit.edu>
#
# Derived from the original crAss v1.1 by:
# (c) Bas E. Dutilh, Robert Schmieder, Jim Nulton, Ben Felts, Peter Salamon, Robert A. Edwards and John L. Mokili
# Please cite our paper "Reference-independent comparative metagenomics using cross-assembly: crAss"
# See crAss web site for more info http://edwards.sdsu.edu/crass/
#
# LICENSE: GPLv3

my %size_options = ();
++$size_options{"minimum"};
++$size_options{"shot"};
++$size_options{"wootters"};
++$size_options{"reads"};

my $usage = "$0 working_directory\n";
my $crAss_dir = $0;
$crAss_dir =~ s/crAss\.pl//;
my $delete_tmp_files = 1;             # 1=YES ; 0=NO
my $print_warnings = 0;               # 1=YES ; 0=NO

# Check input
if (scalar @ARGV != 1) {
    die $usage; }
if (! (-d $ARGV[0])) {
    die "ERROR\t$ARGV[0] is not a directory\n$usage"; }

# Read names for samples
my %names = ();
my %used_name = ();
if (open (IN, "<$ARGV[0]/names.txt")) {
    while (my $line = <IN>) {
	chomp ($line);
	my @a = split /\t/o, $line;
	$a[1] =~ s/\W+/./og;
	if (exists ($used_name{$a[1]})) {
	    die "ERROR\t$a[1] is not a unique name (used for $used_name{$a[1]} and $a[0])\n"; }
	$used_name{$a[1]} = $a[0];
	$names{$a[0]} = $a[1]; } }

# Read files from input directory
opendir (DIR, "$ARGV[0]");
my @files = readdir (DIR);
closedir (DIR);
my @reads_files = ();
my $ace_file = "";
FILE: foreach my $file (@files) {
    if ($file =~ /\.ace\Z/o) {
	if ($ace_file ne "") {
	    die "ERROR\tTwo .ace files in directory $ARGV[0]: $ace_file and $file\n"; }
	$ace_file = $file; }
    elsif (($file =~ /\.fa\Z/i)
	   || ($file =~ /\.fna\Z/i)
	   || ($file =~ /\.fasta\Z/i)
	   || ($file =~ /\.fq\Z/i)
	   || ($file =~ /\.fastq\Z/i)) {
	if ((scalar (keys %names) == 0) || (exists ($names{$file}))) {
	    push (@reads_files, $file); } } }
if (scalar (keys %names) == 0) {
    foreach my $file (@reads_files) {
	my $name = $file;
	$name =~ s/\W+/./og;
	if (exists ($used_name{$name})) {
	    die "ERROR\t$name is not a unique name (used for $used_name{$name} and $file)\n"; }
	$used_name{$name} = $file;
	$names{$file} = $name; } }
my $n_samples = scalar @reads_files;
if ($n_samples < 2) {
    die "ERROR\tNot enough reads files found in $ARGV[0] (only $n_samples, should be at least 2)\n"; }

@reads_files = sort { $names{$a} cmp $names{$b}; } @reads_files;
# RPKM: Store the lengths of the reads.
my %read_lengths = ();

# Read read IDs
my %read2sample = ();
for (my $f = 0; $f < scalar @reads_files; ++$f) {
    my $total_length = 0;
    if (! (open (IN, "<$ARGV[0]/$reads_files[$f]"))) {
	die "ERROR\tCan't open $ARGV[0]/$reads_files[$f]: $!\n"; }
    while (my $line = <IN>) {
	if ($line =~ /^>(\S+)/o) {
	    if (exists ($read2sample{$1})) {
		die "ERROR\tRead ID $1 present in 2 files: $read2sample{$1} and $reads_files[$f]\n"; }
	    $read2sample{$1} = $f; }
	elsif ($line =~ /^\@(\S+)/o) {
	    if (exists ($read2sample{$1})) {
		die "ERROR\tRead ID $1 present in 2 files: $read2sample{$1} and $reads_files[$f]\n"; }
	    $read2sample{$1} = $f; }
	else {
	    # semenko: added length calculations
	    # Validated with: grep -Fv ">" | tr -d '\n' | wc -m
	    chomp($line);
	    $total_length += length($line);
	}
    }
    # RPKM: We downscale these reads by 1 million (the M in RPKM)
    $read_lengths{$f} = $total_length / 1000000;
}
close (IN);

# Read contigs file
if (! (open (IN, "<$ARGV[0]/$ace_file"))) {
    die "ERROR\tCan't open $ARGV[0]/$ace_file: $!\n"; }
my %contig2samples = ();
my %sample2contigs = ();
my %assembled_reads = ();
my %contig_lengths = ();
my $this_contig = "";
while (my $line = <IN>) {
    if ($line =~ /^CO (\S+) (\d+)/o) {
	$this_contig = $1;
	# We downscale these contig lengths by 1000 (the K in RPKM)
	$contig_lengths{$this_contig} = $2 / 1000; }
    elsif ($line =~ /^AF (.+?)(\.\d+\-\d+(\.fm\d+)?(\.to\d+)?)? /o) {
	++$assembled_reads{$1};
	if (! (exists ($read2sample{$1}))) {
	    if ($print_warnings) {
		print STDERR "WARNING\tRead $1 in $ace_file but not in any of the reads files in directory $ARGV[0] - skipping\n"; } }
	else {
	    ++$contig2samples{$this_contig}{$read2sample{$1}};
	    ++$sample2contigs{$read2sample{$1}}{$this_contig}; } } }
close (IN);

# Output reads per contig table
if (! (open (OUT, ">$ARGV[0]/output.contigs2reads.txt"))) {
    die "ERROR\tCan't write to $ARGV[0]/output.contigs2reads.txt: $!\n"; }
for (my $i = 0; $i < scalar @reads_files; ++$i) {
    if (! (exists ($names{$reads_files[$i]}))) {
	die "ERROR\tNo name for sample $reads_files[$i]\n"; }
    print OUT "\t$names{$reads_files[$i]}"; }
print OUT "\n";
my %max = ();
my %n_reads_in_contigs = ();


# Note: This has been modified to support RPKM normalization of contig reads/hits.  -- semenko

# For RPKM, we are going to normalize by read length (per million) and contig length (per thousand).
# Namely, we need to get the # of nucleotides in the .fa files (divided by 1 mil), and the #s of
# nucleotides in each contig in the ACE file (divided by 1,000) and then down-scale the hits.

foreach my $sample (0..scalar(@reads_files)-1) {
    $max{$sample} = 0; }
foreach my $contig (sort keys %contig2samples) {
    my @values = ();
    foreach my $sample (0..scalar(@reads_files)-1) {
	if (! (exists ($contig2samples{$contig}{$sample}))) {
	    $contig2samples{$contig}{$sample} = 0; }

	# First, let's massively scale up the hit value to avoid floating point issues.
	# We've pre-scaled 
	my $normalized_counts = $contig2samples{$contig}{$sample};
	# print "The contig: $contig has length: $contig_lengths{$contig}\n";
	# print "The sample: $sample has length: $read_lengths{$sample}\n";
	# We divide by contig_lengths & read_lengths
	# ** NOTE: Both of these have been already downscaled above. **
	$normalized_counts = $normalized_counts / $contig_lengths{$contig} / $read_lengths{$sample};
	# Optionally, you could take a log to make the heatmap prettier
	$normalized_counts = log($normalized_counts+0.00000000001) / log(2);
	# print "Sample: $sample, $normalized_counts\n";
	push (@values, $normalized_counts);
	# push (@values, $contig2samples{$contig}{$sample});

	# This is used for the "wootters" distance (?)
	$n_reads_in_contigs{$sample} += $contig2samples{$contig}{$sample};
	# Compute a maximum value, used later.
	if ($max{$sample} < $contig2samples{$contig}{$sample}) {
	    $max{$sample} = $contig2samples{$contig}{$sample}; } }
    print OUT "$contig\t", (join ("\t", @values)), "\n"; }

# At the end of the output.contigs2reads.txt file we also report the unassembled reads
# RPKM (@semenko): Commented this out, since it seems useless.
#foreach my $read (keys %read2sample) {
#    if (! (exists ($assembled_reads{$read}))) {
#	print OUT "$read";
#	foreach my $sample (0..scalar(@reads_files)-1) {
#	    if ($sample eq $read2sample{$read}) {
#		print OUT "\t1"; }
#	    else {
#		print OUT "\t0"; } }
#	print OUT "\n"; } }
close (OUT);

# Rounded max number of reads for each sample calculated as max axis value for the 2D/3D plots
foreach my $sample (0..scalar(@reads_files)-1) {
    $max{$sample} = 100 * (int ($max{$sample} / 100) + 1); }

# Output the distance matrix
foreach my $size_correction (keys %size_options) {
    if (! (open (OUT, ">$ARGV[0]/output.$size_correction.distance_matrix.txt"))) {
	die "ERROR\tCan't write to $ARGV[0]/output.$size_correction.distance_matrix.txt: $!\n"; }
    print OUT "  $n_samples\n";
    my $sqrt2 = sqrt (2);
    foreach my $sample_i (0..scalar(@reads_files)-1) {
	print OUT "$names{$reads_files[$sample_i]}";
	my $contigs_i = 0;
	my $reads_i = 0;
	foreach my $contig_i (keys %{$sample2contigs{$sample_i}}) {
	    ++$contigs_i;
	    $reads_i += $contig2samples{$contig_i}{$sample_i}; }
	if ($contigs_i == 0) {
	    die "ERROR\tNone of the contigs contain any reads from $reads_files[$sample_i]\n"; }
	foreach my $sample_j (0..scalar(@reads_files)-1) {
	    my $contigs_j = 0;
	    my $reads_j = 0;
	    foreach my $contig_j (keys %{$sample2contigs{$sample_j}}) {
		++$contigs_j;
		$reads_j += $contig2samples{$contig_j}{$sample_j}; }
	    if ($contigs_j == 0) {
		die "ERROR\tNone of the contigs contain any reads from $reads_files[$sample_j]\n"; }
	    my $contigs_ij = 0;
	    my $reads_ij = 0;
	    my $reads_ji = 0;
	    foreach my $contig_i (keys %{$sample2contigs{$sample_i}}) {
		if (exists ($sample2contigs{$sample_j}{$contig_i})) {
		    ++$contigs_ij;
		    $reads_ij += $contig2samples{$contig_i}{$sample_i};
		    $reads_ji += $contig2samples{$contig_i}{$sample_j}; } }
# minimum-corrected distance
	    my $dist = 1 - ($contigs_ij / (&min ($contigs_i, $contigs_j)));
	    if ($size_correction eq "shot") {
# shot-corrected distance
		$dist = 1 - ($contigs_ij /
			     (($sqrt2 * $contigs_i * $contigs_j) / (sqrt ($contigs_i ** 2 + $contigs_j ** 2)))); }
	    elsif ($size_correction eq "wootters") {
# wootters distance
		$dist = 0;
		foreach my $contig (keys %contig2samples) {
		    $dist += sqrt (($contig2samples{$contig}{$sample_i} / $n_reads_in_contigs{$sample_i}) * ($contig2samples{$contig}{$sample_j} / $n_reads_in_contigs{$sample_j})); }
		if (($dist > 1) || ($sample_i eq $sample_j)) {
		    $dist = 1; }
		$dist = acos ($dist) / (pi / 2); }
	    elsif ($size_correction eq "reads") {
# percentage of reads in shared contigs
		$dist = 1 - sqrt ( .5 * ((($reads_ij ** 2) / ($reads_i ** 2)) + (($reads_ji ** 2) / ($reads_j ** 2)))); }
	    elsif ($size_correction ne "minimum") {
		die "ERROR\tIllegal size correction option: $size_correction\n"; }
	    print OUT "\t$dist"; }
	print OUT "\n"; }
    close (OUT); }

# Output image
#get labels
my (@labels,$tmplabel);
foreach my $n (@reads_files) {
    $tmplabel = $names{$n};
    $tmplabel = substr ($tmplabel, 0, 50).'...' if (length ($tmplabel) > 50);
    $tmplabel =~ s/\_/-/g;
    push (@labels, $tmplabel); }

my $min = 0.9; # minimum value on axes for log plot

# 2D plot for 2 samples
if ($n_samples == 2) {
    if (! (open (OUT, ">$ARGV[0]/tmp.2d.data"))) {
	die "ERROR\tCan't write to $ARGV[0]/tmp.2d.data: $!\n"; }
    print OUT "$names{$reads_files[0]}";
    foreach my $sample (@reads_files) {
	print OUT "\t$names{$sample}"; }
    print OUT "\n";

    foreach my $contig (sort keys %contig2samples) {
	my @values = ();
	foreach my $sample (0..scalar(@reads_files)-1) {
	    if ($contig2samples{$contig}{$sample} == 0) {
		push (@values, $min); }
	    else {
		push (@values, $contig2samples{$contig}{$sample}); } }
	print OUT (join ("\t", @values)), "\n"; }
    close (OUT);
    my $all_max = $min;
    foreach my $sample (0..scalar(@reads_files)-1) {
	if ($all_max < $max{$sample}) {
	    $all_max = $max{$sample}; } }
    if (! (open (GNU, ">$ARGV[0]/tmp.gnuplot.txt"))) {
	die "ERROR\tCan't write to $ARGV[0]/tmp.gnuplot.txt: $!\n"; }
    print GNU "set terminal png nocrop enhanced size 700,660\n";
    print GNU "unset key\n";
    print GNU "set output \"$ARGV[0]/output.image.png\"\n";
    print GNU "set ticslevel 0\n";
    print GNU "set log xy\n";
    print GNU "set xrange [$min:$all_max]\n";
    print GNU "set yrange [$min:$all_max]\n";
    print GNU "set xlabel \"".shift(@labels)."\"\n";
    print GNU "set ylabel \"".shift(@labels)."\"\n";
    print GNU "plot '$ARGV[0]/tmp.2d.data' using 1:2 with points pt 13 lc rgb \"#56B4E9\" notitle\n";
    #print GNU "plot \"$ARGV[0]/tmp.2d.data\" using 1:2 with points notitle\n";
    close (GNU);
    system("gnuplot < $ARGV[0]/tmp.gnuplot.txt"); }

# 3D plot for 3 samples
elsif ($n_samples == 3) {
    if (! (open (OUT, ">$ARGV[0]/tmp.no.projection"))) {
	die "ERROR\tCan't write to $ARGV[0]/tmp.no.projection: $!\n"; }
#	print OUT (join ("\t", @reads_files)), "\n";
    foreach my $contig (sort keys %contig2samples) {
	my @values = ();
	foreach my $sample (0..scalar(@reads_files)-1) {
	    push (@values, $contig2samples{$contig}{$sample}); }
	print OUT (join ("\t", @values)), "\n"; }
    close (OUT);
    foreach my $sample (0..scalar(@reads_files)-1) {
	if (! (open (OUT, ">$ARGV[0]/tmp.$names{$reads_files[$sample]}.projection"))) {
	    die "ERROR\tCan't write to $ARGV[0]/tmp.$names{$reads_files[$sample]}.projection: $!\n"; }
#		print OUT (join ("\t", @reads_files)), "\n";
	foreach my $contig (sort keys %contig2samples) {
	    my @values = ();
	    foreach my $sample_2 (0..scalar(@reads_files)-1) {
		if ($sample eq $sample_2) {
		    push (@values, $min); }
		else {
		    push (@values, $contig2samples{$contig}{$sample_2}); } }
	    print OUT (join ("\t", @values)), "\n"; }
	close (OUT); }
    my $all_max = $min;
    foreach my $sample (0..scalar(@reads_files)-1) {
	if ($all_max < $max{$sample}) {
	    $all_max = $max{$sample}; } }
    if (! (open (GNU, ">$ARGV[0]/tmp.gnuplot.txt"))) {
	die "ERROR\tCan't write to $ARGV[0]/tmp.gnuplot.txt: $!\n"; }
    print GNU "set terminal png nocrop enhanced size 700,660\n";
    #print GNU "unset key\n";
    print GNU "set output \"$ARGV[0]/output.image.png\"\n";
    print GNU "set view 70,45\n";
    print GNU "set ticslevel 0\n";
    print GNU "set log xyz\n";
    print GNU "set xrange [$min:$all_max]\n";
    print GNU "set yrange [$min:$all_max] reverse\n";
    print GNU "set zrange [$min:$all_max]\n";
    print GNU "set xlabel \"".$labels[0]."\"\n";
    print GNU "set ylabel \"".$labels[1]."\"\n";
    print GNU "set zlabel \"".$labels[2]."\"\n";
    print GNU "set border 127+256+512\n";
    print GNU "splot '$ARGV[0]/tmp.no.projection' with points pt 13 lc rgb \"#56B4E9\" notitle, ";
    #print GNU "splot \"$ARGV[0]/tmp.no.projection\" with points notitle, ";
    print GNU "\"$ARGV[0]/tmp.$names{$reads_files[0]}.projection\" with dots lc rgb \"blue\" title \"$labels[0]\", ";
    print GNU "\"$ARGV[0]/tmp.$names{$reads_files[1]}.projection\" with dots lc rgb \"blue\" title \"$labels[1]\", ";
    print GNU "\"$ARGV[0]/tmp.$names{$reads_files[2]}.projection\" with dots lc rgb \"blue\" title \"$labels[2]\"\n";
    close (GNU);
    system ("gnuplot < $ARGV[0]/tmp.gnuplot.txt"); }

# Cladogram for >3 samples
else {
    foreach my $size_correction (keys %size_options) {
# Call BioNJ to generate cladogram
	if (! (open (OUT, ">$ARGV[0]/tmp.BioNJ.pars"))) {
	    die "ERROR\tCan't write to $ARGV[0]/tmp.BioNJ.pars: $!\n"; }
	print OUT "$ARGV[0]/output.$size_correction.distance_matrix.txt\n$ARGV[0]/output.$size_correction.cladogram.ph";
	close (OUT);
	system (($crAss_dir ? $crAss_dir : '.')."/BIONJ_linux < $ARGV[0]/tmp.BioNJ.pars 2>>$ARGV[0]/crass.log");

# Create tree with Drawtree from Phylip
# V - output device
# L - Label align
# D - don't overlap labels
# M - Margins
# F - Font type
# C - font height (relative)
# Y - Finish
	if(-e "$ARGV[0]/output.$size_correction.cladogram.ph" && 0) {
	    #use branch length
	    system ("cd $ARGV[0];rm -f plotfile;echo \"output.$size_correction.cladogram.ph\\nV\\nN\\nL\\nA\\nD\\nM\\n0\\n0\\nF\\nHelvetica\\nC\\n0.5\\nY\\n\" | phylip drawtree") == 0
		or die "ERROR: phylip drawtree (with branch length) system call failed: $?";
	    #create png file from tree using ghostscript
	    if (-e "$ARGV[0]/plotfile") {
		system ("cd $ARGV[0];gs -q -dBATCH -dNOPAUSE -r72 -sOutputFile=output.$size_correction.image.png -sDEVICE=pnggray -dTextAlphaBits=4 -dGraphicsAlphaBits=4 plotfile") == 0
		    or die "ERROR: ghostscript (with branch length) system call failed: $?"; }
	    else {
		die "ERROR: missing plotfile.\n"; }
	    #ignore branch length
	    system ("cd $ARGV[0];rm -f plotfile;echo \"output.$size_correction.cladogram.ph\\nV\\nN\\nL\\nA\\nD\\nM\\n0\\n0\\nF\\nHelvetica\\nC\\n0.5\\nB\\nY\\n\" | phylip drawtree") == 0
		or die "ERROR: phylip drawtree (no branch length) system call failed: $?";
	    #create png file from tree using ghostscript
	    if (-e "$ARGV[0]/plotfile") {
		system ("cd $ARGV[0];gs -q -dBATCH -dNOPAUSE -r72 -sOutputFile=output.$size_correction.image2.png -sDEVICE=pnggray -dTextAlphaBits=4 -dGraphicsAlphaBits=4 plotfile") == 0
		    or die "ERROR: ghostscript (no branch length) system call failed: $?"; }
	    else {
		die "ERROR: missing plotfile.\n"; } }
	else {
	    # die "ERROR: missing tree data file.\n";
	}
    }
}

# Remove temporary files
#if ($delete_tmp_files) {
#    system ("rm $ARGV[0]/tmp.*");
#    system ("rm $ARGV[0]/plotfile") if(-e "$ARGV[0]/plotfile"); }

print STDERR "done\n";

# Subroutines
sub min () {
    my $min = $_[0];
    foreach my $i (@_) {
	if ($i < $min) {
	    $min = $i; } }
    return ($min); }
