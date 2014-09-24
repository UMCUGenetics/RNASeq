#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(tmpnam);
use Getopt::Long;
use File::Basename;

#this script assumes a directory structure like this:
#+<rundir>
#	+<sample>
#		bam

#run by typing:
#perl read_counts.pl /path/to/rundir /path/to/reference_gtf

sub usage{
    warn <<END;
    
    Usage:
    
    Run by typing:	perl read_counts.pl -input /path/to/rundir -gtf /path/to/reference/gtf
    
    Required params:
    -input		[s]	Input directory [/path/to/rundir]
    -gtf		[s]	GTF file of reference genome [/path/to/reference/gtf]

END
    exit;
}


## ======================================
## Get options
## ======================================

my %opt;
%opt = (
    'help'			=> undef,
    'input'			=> undef,
    'gtf'			=> undef,
    'sambamba_path'		=> '/hpc/local/CentOS6/cog_bioinf/bin/sambamba_v0.4.5',
    'python_path'		=> '/hpc/local/CentOS6/cog_bioinf/Python-2.7.6/bin/python2.7',
    'merge_normalize_script'	=> '/hpc/cog_bioinf/data/annelies/RNA_Seq/scripts/merge_normalize_count_tables_test.r'
);

die usage() if @ARGV == 0;
GetOptions (
    'h|help'		=> \$opt{help},
    'i|input=s@'	=> \$opt{input},
    'gtf=s'		=> \$opt{gtf}
) or die usage();

die usage() if $opt{help};
die usage() unless $opt{input};
die usage() unless $opt{output};

# my $gtf = '/hpc/cog_bioinf/GENOMES/CanFam3.1_full/Canis_familiaris.CanFam3.1.75.gtf';
# my $directory = '/hpc/cog_bioinf/data/PROJECTS/Cf_RNAseq_FrankRiemers';

my $directory = $opt{input};
my $gtf = $opt{gtf};

opendir(DIR, $directory);
my @dirs = grep { !/^\.\.?$/ } readdir DIR;
closedir(DIR);


my $mainJobID = "$directory/QSUB_".get_job_id().".sh";
open QSUB, ">$mainJobID";
print QSUB "\#!/bin/sh\n\#\$ \-o $directory\n\#\$ \-e $directory\n\n";

my @hold_ids=();

foreach my $dir (@dirs){ 
    print "$dir\n";
    my ($HTSeqCount_job_id) = ("htseqCount_".get_job_id());
    push (@hold_ids, $HTSeqCount_job_id);
    open HT,">$directory/$dir/$HTSeqCount_job_id.sh";
    print HT "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    
    my $bam = "$directory/$dir/merged_$dir\_20144007.bam";
    my $sortedBam = "$directory/$dir/merged_$dir\_20144007_sortedByName.bam";
    
    print HT "$opt{sambamba_path} sort -n -t 4 -o $sortedBam $bam\n";
    print HT "$opt{sambamba_path} view $sortedBam | $opt{python_path} -m HTSeq.scripts.count -m union -s yes -i gene_id - $gtf > $directory/$dir/$dir\_read_counts.txt\n";
    
    print HT "rm $sortedBam\n";

    close HT;
    
    print QSUB "qsub -q veryshort -pe threaded 4 -o $directory/$dir -e $directory/$dir -N $HTSeqCount_job_id $directory/$dir/$HTSeqCount_job_id.sh\n\n";
}

my ($normTables_job_id) = ("normTables_".get_job_id());
open MT,">$directory/$normTables_job_id.sh";
print MT "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
print MT "R --vanilla --args $directory < $opt{merge_normalize_script}\n";
close MT;
my $hold_line = join ',',@hold_ids;
print QSUB "qsub -q veryshort -o $directory -e $directory -N $normTables_job_id -hold_jid $hold_line $directory/$normTables_job_id.sh\n\n";

close QSUB;
system "sh $mainJobID";

sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}

