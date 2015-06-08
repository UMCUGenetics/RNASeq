#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(tmpnam);
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);

## ======================================
## Usage: see -h
## ======================================

sub usage{
  warn <<END;
  
  Usage:
  
  Run by typing:	perl RNAseqAnalyse.pl -input [run directory] -outputDir [output directory]
  
    Required params:
    -i|input                          [s]	Input run directory [/path/to/rundir] (list of input run directories is allowed, separated by ',')
    -o|outputDir                      [s]	Directory to store mapping results [/path/to/outputDir]
    
    Options:
    -h|help                           [s]	Help
    -sp|species                       [s]	Species of RNA-seq data [HUMAN/RAT/MOUSE/ZEBRAFISH/DOG] Default: HUMAN
    -pe|nthreads                      [i]	Number of threads for processes run on execute nodes. [number] Default: 4
    -g|genome                         [s]	Directory with reference genome [/path/to/genome] Default: /hpc/cog_bioinf/GENOMES
    -fqc|fastqc                       [s]	Perform FastQC? [yes/no] Default: yes
    -m|mapping                        [s]	Perform mapping with STAR? [yes/no] Default: yes
    -stranded                         [s]	Is the RNA-seq data from a strand-specific assay? [yes/no/reversed] Default: reversed
    -fusionSearch                     [s]	Want to detect fusion genes? [yes/no] Default: yes
    -c|count                          [s]	Want to count the mapped reads? [yes/no] Default: yes
    -n|normalize                      [s]	Want to merge and normalize counted reads? [yes/no] Default: yes
    -rpkm                             [s]	Want to retrieve RPKMs? [yes/no] Default: yes
    -bamqc                            [s]	Want to retrieve bam statistics? [yes/no] Default: yes
    -id                               [s]	Which id attribute do you want to count? [gene_id/transcript_id/exon_id] Default: gene_id
    -uniq                             [s]	Want to retreive alignment file (BAM) with unique reads? [yes/no] Default: no
    -chimSegmentMin                   [i]	Minimum length of chimeric segment length. [number] Default: 15
    -chimJunctionOverhangMin          [i]	Minimum overhang for a chimeric junction. [number] Default: 15
    -outSJfilterIntronMaxVsReadN      [s]	Maximum gap allowed for junctions supported by 1,2,3...N reads. [number] Default: 10.000.000
      
END
  exit;
}



## ======================================
## Get options
## ======================================

my %opt;
%opt = (
    'help'				=> undef,
    'input'				=> undef,
    'outputDir'				=> undef,
    'nthreads'				=> 4,
    'fastqc'				=> "yes",
    'mapping'				=> "yes",
    'stranded'				=> "reversed",
    'fusionSearch'			=> "yes",
    'count'				=> "yes",
    'normalize'				=> "yes",
    'rpkm'				=> "yes",
    'bamqc'				=> "yes",
    'id'				=> "gene_id",
    'uniq'				=> "no",
    'chimSegmentMin'			=> 15,
    'outSJfilterIntronMaxVsReadN'	=> 10000000,
    'chimJunctionOverhangMin'		=> 15,
    'species'				=> "HUMAN",
    'genome'				=> '/hpc/cog_bioinf/GENOMES',
    'fasta'				=> undef,
    'intervallist'			=> undef,
    'gtf_file'				=> undef,
    'refflat_file'			=> undef,
    'genesizes_file'			=> undef,
    'fastqc_path'			=> '/hpc/cog_bioinf/common_scripts/FastQC/fastqc',
#     'star_path'				=> '/hpc/local/CentOS6/cog_bioinf/STAR_2.3.0e/STAR',
    'star_path'				=> '/hpc/local/CentOS6/cog_bioinf/STAR-STAR_2.4.0i/source/STAR',
    'sambamba_path'			=> '/hpc/local/CentOS6/cog_bioinf/bin/sambamba_v0.4.5',
    'picard_path'			=> '/hpc/cog_bioinf/common_scripts/picard-tools-1.98',
    'python_path'			=> '/hpc/local/CentOS6/cog_bioinf/Python-2.7.6/bin/python2.7',
    'bamstats_path'			=> '/hpc/cog_bioinf/common_scripts/bamMetrics/bamMetrics.pl'
);

die usage() if @ARGV == 0;
GetOptions (
    'h|help'				=> \$opt{help},
    'i|input=s@'			=> \$opt{input},
    'o|outputDir=s'			=> \$opt{outputDir},
    'pe|nthreads=i'			=> \$opt{nthreads},
    'g|genome=s'			=> \$opt{genome},
    'fa|fasta=s'			=> \$opt{fasta},
    'stranded=s'			=> \$opt{stranded},
    'fqc|fastqc=s'			=> \$opt{fastqc},
    'm|mapping=s'			=> \$opt{mapping},
    'fusionSearch=s'			=> \$opt{fusionSearch},
    'c|count=s'				=> \$opt{count},
    'n|normalize=s'			=> \$opt{normalize},
    'rpkm=s'				=> \$opt{rpkm},
    'bamqc=s'				=> \$opt{bamqc},
    'id=s'				=> \$opt{id},
    'uniq=s'				=> \$opt{uniq},
    'chimSegmentMin=i'			=> \$opt{chimSegmentMin},
    'outSJfilterIntronMaxVsReadN=s'	=> \$opt{outSJfilterIntronMaxVsReadN},
    'chimJunctionOverhangMin=i'		=> \$opt{chimJunctionOverhangMin},
    'sp|species=s'			=> \$opt{species}
) or die usage();

#check input paramaters
die usage() if $opt{help};
die "[ERROR] Number of threads must be at least 3!\n" if ($opt{nthreads} < 3);
die usage() unless ( $opt{input} );
die usage() unless ( $opt{outputDir} );
die "[ERROR] Nothing to do!\n" if ( ($opt{fastqc} eq "no") && ($opt{mapping} eq "no") && ($opt{count} eq "no") && ($opt{normalize} eq "no") && ($opt{rpkm} eq "no") && ($opt{bamqc} eq "no") );
die "[ERROR] Wrong option ($opt{fastqc}) given with -fastqc. Must be \"yes\" or \"no\".\n" if ( ($opt{fastqc} ne "no") && ($opt{fastqc} ne "yes") );
die "[ERROR] Wrong option ($opt{mapping}) given with -mapping. Must be \"yes\" or \"no\".\n" if ( ($opt{mapping} ne "no") && ($opt{mapping} ne "yes") );
die "[ERROR] Wrong option ($opt{stranded}) given with -stranded. Must be \"yes\", \"no\" or \"reversed\".\n" if ( ($opt{stranded} ne "no") && ($opt{stranded} ne "yes") && ($opt{stranded} ne "reversed") );
die "[ERROR] Wrong option ($opt{fusionSearch}) given with -fusionSearch. Must be \"yes\" or \"no\".\n" if ( ($opt{fusionSearch} ne "no") && ($opt{fusionSearch} ne "yes") );
die "[ERROR] Value given with -chimSegmentMin ($opt{chimSegmentMin}) must be higher than 0, when searching for fusion genes.\n" if ( ($opt{fusionSearch} eq "yes") && ($opt{chimSegmentMin} <= 0) );
die "[ERROR] Wrong option ($opt{count}) given with -count. Must be \"yes\" or \"no\".\n" if ( ($opt{count} ne "no") && ($opt{count} ne "yes") );
die "[ERROR] Wrong option ($opt{normalize}) given with -normalize. Must be \"yes\" or \"no\".\n" if ( ($opt{normalize} ne "no") && ($opt{normalize} ne "yes") );
die "[ERROR] Wrong option ($opt{rpkm}) given with -rpkm. Must be \"yes\" or \"no\".\n" if ( ($opt{rpkm} ne "no") && ($opt{rpkm} ne "yes") );
die "[ERROR] Wrong option ($opt{uniq}) given with -uniq. Must be \"yes\" or \"no\".\n" if ( ($opt{uniq} ne "no") && ($opt{uniq} ne "yes") );
die "[ERROR] Wrong option ($opt{id}) given with -id. Must be \"gene_id\", \"transcript_id\" or \"exon_id\".\n" if ( ($opt{id} ne "gene_id") && ($opt{id} ne "transcript_id") && ($opt{id} ne "exon_id") );

my $SPECIES = uc $opt{species};
if ($SPECIES eq "HUMAN"){
    $opt{genome} .= '/Homo_sapiens.GRCh37.GATK.illumina/STAR';
    $opt{fasta} = $opt{genome}.'/hg19.fa';
    $opt{gtf_file} = $opt{genome}.'/Homo_sapiens.GRCh37.74.gtf';
    $opt{refflat_file} = $opt{genome}.'/hg19.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/Homo_sapiens.GRCh37.GATK.illumina.rRNA.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Homo_sapiens.GRCh37.74_exon_gene_sizes.txt';
} elsif ($SPECIES eq "RAT"){
    $opt{genome} .= '/rat_GATK_illumina_rnor_50/STAR';
    $opt{fasta} = $opt{genome}.'/rat_GATK_illumina_rnor_50/STAR/rnor50.fa';
    $opt{gtf_file} = $opt{genome}.'/Rattus_norvegicus.Rnor_5.0.71.gtf';
    $opt{refflat_file} = $opt{genome}.'/rnor50.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/rno5_rRNA_intervallist.txt';
    $opt{genesizes_file} = $opt{genome}.'/Rattus_norvegicus.Rnor_5.0.71_exon_gene_sizes.txt';
} elsif ($SPECIES eq "MOUSE"){
    $opt{genome} .= '/Mus_musculus_GRCm38_GATK_illumina_bwa075/STAR/mus_musculus_GRCm38';
    $opt{fasta} = $opt{genome}.'/Mm_GRCm38_gatk_sorted.fa';
    $opt{gtf_file} = $opt{genome}.'/Mus_musculus.GRCm38.70.gtf';
    $opt{refflat_file} = $opt{genome}.'/Mus_musculus_GRCm38.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/Mus_musculus_GRCm38.rRNA.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Mus_musculus.GRCm38.70_exon_gene_sizes.txt';
} elsif ($SPECIES eq "ZEBRAFISH"){
    $opt{genome} .= '/zfish9/STAR';
    $opt{fasta} = $opt{genome}.'/Zv9_66.fa';
    $opt{gtf_file} = $opt{genome}.'/Danio_rerio.Zv9.75.gtf';
    $opt{refflat_file} = $opt{genome}.'/zfish9.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/Danio_rerio.Zv9.75.rRNA.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Danio_rerio.Zv9.75_exon_gene_sizes.txt';
} elsif ($SPECIES eq "DOG"){
    $opt{genome} .= '/CanFam3.1.GATK.illumina/STAR';
    $opt{fasta} = $opt{genome}.'/cf3_ens71_GATK.fa';
    $opt{gtf_file} = $opt{genome}.'/Canis_familiaris.CanFam3.1.75.gtf';
    $opt{refflat_file} = $opt{genome}.'/CanFam3.1.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/CanFam3.1_rRNA_genes.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Canis_familiaris.CanFam3.1.75_exon_gene_sizes.txt';
} else {
    die "[ERROR] Wrong species ($SPECIES). Only HUMAN, RAT, MOUSE, DOG or ZEBRAFISH genomes are allowed.\n"
}



## ======================================
## Retrieve inputfiles
## ======================================

my @samplefiles;
my $input = $opt{input};
my @input = split(/,/,join(',',@$input));

foreach my $fastqdir (@input){
    print "$fastqdir\n";
    
    if ( ! -e $fastqdir ){ die "$fastqdir does not exist." }
    
    open (FIND, "find $fastqdir -name '*_R1_001.fastq.gz' |");
    while (my $f= <FIND>) {
	chomp $f;
	my $pattern = 'Undetermined';
	push @samplefiles, $f unless $f =~ /$pattern/;
    }
    close FIND;
}

if (scalar @samplefiles ==0) {
    print "Nothing to do!!\n\n";
    usage;
}



## ======================================
## Create directories
## ======================================

my $rundir = "";
my $samples = {};

foreach my $f (@samplefiles){
        
    my $base = basename($f);
    my @parts = split("_",$base);
    my $sample = $parts[0];
    my @input_parts = split("/","$opt{input}");
    
    $rundir=$opt{outputDir};
    push(@{$samples->{$sample}}, $f);
        
    if(! -e "$rundir"){
	make_path($rundir) or die "Couldn't create directory: $rundir\n";
    }
    if(! -e "$rundir/read_counts" && $opt{count} ne 'no'){
	mkdir("$rundir/read_counts") or die "Couldn't create directory: $rundir/read_counts\n";
    }
    if(! -e "$rundir/$sample"){
	mkdir("$rundir/$sample") or die "Couldn't create directory: $rundir/$sample\n";
    }
    if(! -e "$rundir/$sample/fastqc" && $opt{fastqc} ne 'no'){
        mkdir("$rundir/$sample/fastqc") or die "Couldn't create directory: $rundir/$sample/fastqc\n";
    }
    if(! -e "$rundir/$sample/mapping" && $opt{mapping} ne 'no'){
        mkdir("$rundir/$sample/mapping") or die "Couldn't create directory: $rundir/$sample/mapping\n";
    }
    if(! -e "$rundir/$sample/read_counts" && $opt{count} ne 'no'){
	mkdir("$rundir/$sample/read_counts") or die "Couldn't create directory: $rundir/$sample/read_counts\n";
    }
    if(! -e "$rundir/$sample/jobs"){
        mkdir("$rundir/$sample/jobs") or die "Couldn't create directory: $rundir/$sample/jobs\n";
    }
    if(! -e "$rundir/$sample/logs"){
        mkdir("$rundir/$sample/logs") or die "Couldn't create directory: $rundir/$sample/logs\n";
    }
}

## ======================================
## Create jobs
## ======================================

my $mainJobID = "$rundir/".get_job_id()."_qsub.sh";
open QSUB, ">$mainJobID";
print QSUB "\#!/bin/sh\n\#\$ \-o $rundir\n\#\$ \-e $rundir\n\n";

my @hold_ids=();
my @hold_mapping_ids=();
my $paired = 0;

foreach my $sample (keys %{$samples}) {
	print "\nFiles for sample $sample:\n";
	my $job_id = "STAR_$sample\_".get_job_id();
	push @hold_mapping_ids, $job_id;
	#create bash script for STAR submission of this sample
	open STAR_SH,">$rundir/$sample/jobs/$job_id.sh" or die "Couldn't create $rundir/$sample/jobs/$job_id.sh\n";
	print STAR_SH "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	print STAR_SH "uname -n > $rundir/$sample/logs/$sample.host\n";
	print STAR_SH "echo \"mapping pair\t\" `date` >> $rundir/$sample/logs/$sample.host\n";
	print STAR_SH "mkdir -p $rundir/$sample/mapping/tmp/$job_id/\n";
	print STAR_SH "cd $rundir/$sample/mapping/tmp/$job_id/\n\n";
	#Create STAR command
	my $star_command = "$opt{star_path} --genomeDir $opt{genome} --runThreadN $opt{nthreads} --outFileNamePrefix $sample\_ --outReadsUnmapped Fastx --readFilesCommand zcat";
	if ($opt{stranded} eq "no"){
	    $star_command .= " --outSAMstrandField intronMotif";
	}
	if (defined $opt{outSJfilterIntronMaxVsReadN}) {
	    $star_command .= " --outSJfilterIntronMaxVsReadN " . $opt{outSJfilterIntronMaxVsReadN};
	}
	if (defined $opt{chimJunctionOverhangMin}) {
	    $star_command .= " --chimJunctionOverhangMin " . $opt{chimJunctionOverhangMin};
	}
	if ($opt{fusionSearch} eq "yes"){
	    $star_command .= " --chimSegmentMin $opt{chimSegmentMin}";
	}
	elsif ($opt{fusionSearch} eq "no"){
	    $star_command .= " --chimSegmentMin 0";
	}
	
	my ($ID,$PU,$LB,$SM,$PL)='';
	
	#Find Fastq files
	my $R1_fastqs = "";
	my $R2_fastqs = "";
	my $rest_fastqs = "";
	foreach my $fastq (@{$samples->{$sample}}){
	    my $R1 = $fastq;
	    my $R2 = $fastq;
	    $R2 =~ s/R1_001/R2_001/;
	    
	    my $fastqname = basename($R1);
	    my @cols = split("_",$fastqname);
	    $ID = join("_",@cols[0..2]);
	    $PU = $cols[1];
	    $LB = $cols[2];
	    $SM = $cols[0];
	    $PL = 'ILLUMINA';
	    
	    if (-e $R2) {
		$paired = 1;
		$R1_fastqs .= "$R1,";
		$R2_fastqs .= "$R2,";
	    
		print "Pair found:\n";
		print $R1," = R1\n";
		print $R2," = R2\n\n";
	    }else{
		$rest_fastqs .= "$R1,";
		print "Single tag found\n";
		print $R1," = R1\n";
	    }
	    
	    #FASTQC
	    if ($opt{fastqc} eq "yes") {
		my ($FastQC1_job_id, $FastQC2_job_id) = ("FastQC1_".get_job_id(), "FastQC2_".get_job_id());
		
		
		open FASTQC,">$rundir/$sample/jobs/$FastQC1_job_id.sh";
		print FASTQC "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
		print FASTQC "cd $rundir/$sample/fastqc\n\n";
		print FASTQC "uname -n > ../logs/$FastQC1_job_id.host\n";
		print FASTQC "echo \"FastQC\t\" `date` >> ../logs/FastQC_$sample.host\n";
		print FASTQC "$opt{fastqc_path} $R1 -o $rundir/$sample/fastqc\n";
		print QSUB "qsub -q veryshort -o $rundir/$sample/logs -e $rundir/$sample/logs -N $FastQC1_job_id $rundir/$sample/jobs/$FastQC1_job_id.sh\n";
		close FASTQC;
		
		
		if ($paired == 1) {
		    open FASTQC2,">$rundir/$sample/jobs/$FastQC2_job_id.sh";
		    print FASTQC2 "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
		    print FASTQC2 "cd $rundir/$sample/fastqc\n\n";
		    print FASTQC2 "uname -n > ../logs/$FastQC2_job_id.host\n";
		    print FASTQC2 "echo \"FastQC\t\" `date` >> ../logs/FastQC2_$sample.host\n";
		    print FASTQC2 "$opt{fastqc_path} $R2 -o $rundir/$sample/fastqc\n";
		    print QSUB "qsub -q veryshort -o $rundir/$sample/logs -e $rundir/$sample/logs -N $FastQC2_job_id $rundir/$sample/jobs/$FastQC2_job_id.sh\n";
		    close FASTQC2;
		}
	    }
	}
		
	chop $R1_fastqs;
	chop $R2_fastqs;
	chop $rest_fastqs;
	#add fastq files to STAR command
	$star_command .= " --readFilesIn $R1_fastqs $R2_fastqs $rest_fastqs 1>>$rundir/$sample/logs/$sample\_star.log 2>>$rundir/$sample/logs/$sample\_star.err \n\n" unless ( $opt{stranded} eq "reversed" );
	$star_command .= " --readFilesIn $R2_fastqs $R1_fastqs $rest_fastqs 1>>$rundir/$sample/logs/$sample\_star.log 2>>$rundir/$sample/logs/$sample\_star.err \n\n" if ( $opt{stranded} eq "reversed" );
	
	#Print STAR command to star submission script
	print STAR_SH $star_command;

	#Sambamba
	#sam to bam
	print STAR_SH "$opt{sambamba_path} view -t $opt{nthreads} -S -f bam $rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.sam > $rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.bam\n";
	
	#Add read groups using picard
	print STAR_SH "java -jar $opt{picard_path}\/AddOrReplaceReadGroups.jar INPUT=$rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.bam OUTPUT=$rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.rgAdded.bam RGID=$ID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM\n";

	#sort & index bam
	print STAR_SH "$opt{sambamba_path} sort -t $opt{nthreads} -o $rundir/$sample/mapping/$sample\_sorted.bam $rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.rgAdded.bam\n";
	print STAR_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sorted.bam\n";

	#sort bam by name (needed for htseq count)
# 	print STAR_SH "$opt{sambamba_path} sort -n -t $opt{nthreads} -o $rundir/$sample/mapping/$sample\_sorted_byName.bam $rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.rgAdded.bam\n";
	
	if ( $opt{uniq} eq 'yes' ){
	    print STAR_SH "$opt{sambamba_path} view -t $opt{nthreads} -f bam -F \"not secondary_alignment\" $rundir/$sample/mapping/$sample\_sorted.bam >$rundir/$sample/mapping/$sample\_sorted_uniq.bam\n";
	    print STAR_SH "java -Xmx2g  -jar $opt{picard_path}\/MarkDuplicates.jar INPUT=$rundir/$sample/mapping/$sample\_sorted_uniq.bam OUTPUT=$rundir/$sample/mapping/$sample\_sorted_uniq_markDup.bam METRICS_FILE=$rundir/$sample/mapping/$sample\_sorted_uniq_markDup_metrics.txt\n";
	    print STAR_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sorted_uniq.bam\n";
	    print STAR_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sorted_uniq_markDup.bam\n";
	    print STAR_SH "$opt{sambamba_path} flagstat $rundir/$sample/mapping/$sample\_sorted_uniq_markDup.bam\n";
	}

	#remove and move files
	print STAR_SH "rm $rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.sam\n";
	print STAR_SH "rm $rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.bam\n";
	print STAR_SH "rm $rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.rgAdded.bam\n";
	print STAR_SH "mv $rundir/$sample/mapping/tmp/$job_id/* $rundir/$sample/mapping/\n";
	print STAR_SH "rm -r $rundir/$sample/mapping/tmp\n";
	print STAR_SH "rm -r $rundir/$sample/mapping/$sample\__tmp\n";
	
	close STAR_SH;
	
	print QSUB "\n";
	
	if ($opt{mapping} eq "yes") {
	    print QSUB "qsub -q short -pe threaded $opt{nthreads} -o $rundir/$sample/logs -e $rundir/$sample/logs -N $job_id $rundir/$sample/jobs/$job_id.sh\n\n";
	}
	
	#get read counts of mapped reads
	if ($opt{count} eq "yes") {
	    my ($HTSeqCount_job_id) = ("htseqCount_".get_job_id());
	    push (@hold_ids, $HTSeqCount_job_id);
	    open HT,">$rundir/$sample/jobs/$HTSeqCount_job_id.sh";
	    print HT "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	    print HT "uname -n > $rundir/$sample/logs/htseqCount_$sample.host\n";
	    
# 	    my $file = "$rundir/$sample/mapping/$sample\_sorted_byName.bam";
# 	    if ($opt{mapping} eq "no" && ! -e $file) {
# 		print HT "$opt{sambamba_path} sort -n -t $opt{nthreads} -o $rundir/$sample/mapping/$sample\_sorted_byName.bam $rundir/$sample/mapping/$sample\_sorted.bam\n";
# 	    }
	    
	    #set correct strandedness for htseq count
	    my $s_val;
	    if ( $paired==1 ){
		$s_val = 'yes' if $opt{stranded} eq "reversed";
		$s_val = 'reverse' if $opt{stranded} eq "yes";
		$s_val = 'no' if $opt{stranded} eq "no";
	    }elsif ( $paired==0 ){
		$s_val = 'yes' if $opt{stranded} eq "yes";
		$s_val = 'reverse' if $opt{stranded} eq "reversed";
		$s_val = 'no' if $opt{stranded} eq "no";
	    }
	    
# 	    print HT "$opt{sambamba_path} view $rundir/$sample/mapping/$sample\_sorted_byName.bam | $opt{python_path} -m HTSeq.scripts.count -m union -s $s_val -i $opt{id} - $opt{gtf_file} > $rundir/$sample/read_counts/$sample\_htseq_counts.txt\n";
	    print HT "$opt{sambamba_path} view $rundir/$sample/mapping/$sample\_sorted.bam | $opt{python_path} -m HTSeq.scripts.count -m union -r pos -s $s_val -i $opt{id} - $opt{gtf_file} > $rundir/$sample/read_counts/$sample\_htseq_counts.txt\n";

# 	    print HT "rm $rundir/$sample/mapping/$sample\_sorted_byName.bam*\n";
	    close HT;
	    if ( $opt{mapping} eq "no" ){
		print QSUB  "qsub -q veryshort -o $rundir/$sample/logs -e $rundir/$sample/logs -N $HTSeqCount_job_id $rundir/$sample/jobs/$HTSeqCount_job_id.sh\n\n";
	    }else{
		print QSUB  "qsub -q veryshort -o $rundir/$sample/logs -e $rundir/$sample/logs -N $HTSeqCount_job_id -hold_jid $job_id $rundir/$sample/jobs/$HTSeqCount_job_id.sh\n\n";
	    }
	}
}

if ( $opt{normalize} eq "yes"){
    #merge count tables of all samples and normalize the merged table
    my ($normTables_job_id) = ("normTables_".get_job_id());
    open MT,">$rundir/$normTables_job_id.sh";
    print MT "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print MT "uname -n > $rundir/$normTables_job_id.host\n";
    print MT "export MODULEPATH=\$HOME/modules:/hpc/local/CentOS6/cog_bioinf/modules:\${MODULEPATH}\n";
    print MT "module load R_3.1.2\n";
    normalizeRscript($rundir);
    print MT "time R --save < $rundir/normalize.R\n";
    print MT "module unload R_3.1.2\n";
    close MT;
    my $hold_line = join ',',@hold_ids;
    print QSUB "qsub -q veryshort -o $rundir -e $rundir -N $normTables_job_id -hold_jid $hold_line $rundir/$normTables_job_id.sh\n\n";
    push (@hold_ids, $normTables_job_id);
}

if( $opt{rpkm} eq "yes"){
    my ($rpkm_job_id) = ("rpkm_".get_job_id());
    open RP, ">$rundir/$rpkm_job_id.sh";
    print RP "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print RP "uname -n > $rundir/$rpkm_job_id.host\n";
    print RP "export MODULEPATH=\$HOME/modules:/hpc/local/CentOS6/cog_bioinf/modules:\${MODULEPATH}\n";
    print RP "module load R_3.1.2\n";
    rpkmRscript($rundir);
    print RP "time R --save --args $rundir < $rundir/rpkm.R\n";
    print RP "module unload R_3.1.2\n";
    close RP;
    my $hold_line = join ',',@hold_ids;
    print QSUB "qsub -q veryshort -o $rundir -e $rundir -N $rpkm_job_id -hold_jid $hold_line $rundir/$rpkm_job_id.sh\n\n";    
}

if ( $opt{bamqc} eq "yes"){
    my ($bamQC_job_id) = ("bamQC_".get_job_id());
    open BQ, ">$rundir/$bamQC_job_id.sh";
    print BQ "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
    print BQ "uname -n > $rundir/$bamQC_job_id.host\n";

    my $picard_strand;
    if ( $paired==1 ){
	$picard_strand = 'FIRST_READ_TRANSCRIPTION_STRAND' if $opt{stranded} eq 'reversed';
	$picard_strand = 'SECOND_READ_TRANSCRIPTION_STRAND' if $opt{stranded} eq 'yes';
	$picard_strand = 'NONE' if $opt{stranded} eq 'no';
    }elsif ( $paired==0 ){
	$picard_strand = 'FIRST_READ_TRANSCRIPTION_STRAND' if $opt{stranded} eq 'yes';
	$picard_strand = 'SECOND_READ_TRANSCRIPTION_STRAND' if $opt{stranded} eq 'reversed';
	$picard_strand = 'NONE' if $opt{stranded} eq 'no';
    }
    
    my $runname = basename( $rundir );
    
    print BQ "shopt -s nullglob\n\n";
    print BQ "filearray=(\"$rundir\"/*/mapping/*_sorted.bam)\n\n";
    print BQ "bar=\$(printf \",%s\" \"\${filearray[\@]}\")\n\n";
    print BQ "bamline=`echo \$bar | sed 's/,/ -bam /g'`\n\n";
    if ( $paired==0 ){
	print BQ "perl $opt{bamstats_path} \${bamline} -debug -rna -ref_flat $opt{refflat_file} -ribosomal_intervals $opt{intervallist} -strand $picard_strand -single_end -genome $opt{fasta} -run_name $runname -output_dir $rundir/bamMetrics\n\n";
    }
    elsif ( $paired==1 ){
	print BQ "perl $opt{bamstats_path} \${bamline} -debug -rna -ref_flat $opt{refflat_file} -ribosomal_intervals $opt{intervallist} -strand $picard_strand -genome $opt{fasta} -run_name $runname -output_dir $rundir/bamMetrics\n\n";
    }
    close BQ;
    if ( $opt{mapping} eq "no" ){
	print QSUB "qsub -q veryshort -o $rundir -e $rundir -N $bamQC_job_id $rundir/$bamQC_job_id.sh\n\n";
    }else{
	my $hold_line = join ',',@hold_mapping_ids;
	print QSUB "qsub -q veryshort -o $rundir -e $rundir -N $bamQC_job_id -hold_jid $hold_line $rundir/$bamQC_job_id.sh\n\n";
    }
}


close QSUB;
system "sh $mainJobID";



## SUBROUTINES ##

sub normalizeRscript{
    my $rundir = shift;
    my $runname = basename($rundir);
    open normRscript, ">$rundir/normalize.R" or die "cannot open Rscript\n";
    print normRscript <<EOS;
library("DESeq")
dirs <- list.dirs(path="$rundir",recursive=F,full.names=F)
nr_cols=0
sampledir = ''
for ( dir in dirs ){
    sample=basename(dir)
    if ( (sample != 'read_counts') && (sample != 'bamMetrics') ){
	nr_cols = nr_cols+1
	sampledir = dir
    }
}
nr_rows <- nrow(read.table(paste("$rundir/",sampledir,'/read_counts/',basename(sampledir),'_htseq_counts.txt',sep=""),sep="\\t",header=F))-5

output <- matrix(ncol=nr_cols+1, nrow=nr_rows)
col_count=1
samplenames=c('gene')

for ( dir in dirs ) {
    sample=basename(dir)
    if ( (sample != 'read_counts') && (sample != 'bamMetrics') ) {
	samplenames <- append(samplenames, as.character(sample))
	countsfile=paste("$rundir/",sample,'/read_counts/',sample,'_htseq_counts.txt',sep="")
	counts=read.table(countsfile,sep="\\t",header=F)
	colnames(counts) <- c('gene',sample)
	genes<-counts\$gene

	if ( col_count == 1 ){
	    output[,col_count] <- as.character(genes[1:(length(genes)-5)])
	    col_count=col_count+1
	}
	
	output[,col_count] <- as.numeric(counts[1:(nrow(counts)-5),2])
	col_count=col_count+1
    }
}

merged_table <- as.data.frame(output)
colnames(merged_table) <- samplenames

outfile1="$rundir/read_counts/$runname\_readCounts_raw.txt"
write.table(merged_table,file=outfile1,row.names=F,col.names=T,quote=F,sep="\\t")



#normalize raw read counts

counts=read.delim(outfile1, header=T, row.names=1)

conds <- factor(c(colnames(counts)))
cds <- newCountDataSet( counts, conds )
cds <- estimateSizeFactors( cds )

normalized_counts <- round(counts(cds, normalized=TRUE))

outfile2="$rundir/read_counts/$runname\_readCounts_normalized.txt"
write.table(normalized_counts, file=outfile2, row.names=T,col.names=NA, quot=F,sep="\\t")

EOS
    close(normRscript);
}

sub rpkmRscript{
    my $rundir = shift;
    my $runname = basename($rundir);
    my $gene_sizes = $opt{genesizes_file};
    open rpkmRscript, ">$rundir/rpkm.R" or die "cannot open Rscript\n";
    print rpkmRscript <<EOS;
library(edgeR)

exon_gene_sizes <- read.table("$gene_sizes",sep="\\t",header=F,row.names=1)
raw_read_counts <- read.table("$rundir/read_counts/$runname\_readCounts_raw.txt", sep="\\t",header=T,row.names=1)

nrsamples <- ncol(raw_read_counts)
nrrows <- nrow(raw_read_counts)

tab <- matrix(data=NA, nrow=nrrows, ncol=nrsamples)

for (j in 1:nrsamples){
    RPKM = rpkm(raw_read_counts[j], exon_gene_sizes, normalized.lib.sizes=F, log=F)
    for(i in 1:nrrows){
	tab[i,j] = RPKM\$V2[i]
    }
}

df <- data.frame(tab)
colnames(df) <- colnames(raw_read_counts)
rownames(df) <- rownames(raw_read_counts)

outfile="$rundir/read_counts/$runname\_readCounts_RPKM.txt"
write.table(df, file=outfile, row.names=T,col.names=NA, quot=F,sep="\t")

EOS
    close(rpkmRscript);
}

sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}

