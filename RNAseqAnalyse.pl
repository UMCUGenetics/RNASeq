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
  
  Run by typing:	perl RNAseqAnalyse.pl -input [run directory] -output [output directory]
  
    Required params:
    -i|input				[s]	Input run (!) directory [/path/to/rundir]
    -o|output				[s]	Directory to store mapping results [/path/to/store]
    
    Options:
    -h|help				[s]	Help
    -sp|species				[s]	Species of RNA-seq data [HUMAN/RAT/MOUSE/ZEBRAFISH] Default: HUMAN
    -pe|nthreads			[i]	Number of threads for processes run on execute nodes. [number] Default: 4
    -g|genome				[s]	Directory with reference genome [/path/to/genome] Default: /hpc/cog_bioinf/GENOMES
    -fqc|fastqc				[s]	Perform FastQC? [yes/no] Default: yes
    -m|mapping				[s]	Perform mapping with STAR? [yes/no] Default: yes
    -stranded				[s]	Is the RNA-seq data from a strand-specific assay? [yes/no/reversed] Default: reversed
    -fusionSearch			[s]	Want to detect fusion genes? [yes/no] Default: yes
    -c|count				[s]	Want to count the mapped reads? [yes/no] Default: yes
    -feature				[s]	Which feature do you want to count? [gene_id/transcript_id/exon_id] Default: gene_id
    -uniq				[s]	Want to retreive alignment file (BAM) with unique reads? [yes/no] Default: no
    -chimSegmentMin			[i]	Minimum length of chimeric segment length. [number] Default: 15
    -chimJunctionOverhangMin		[i]	Minimum overhang for a chimeric junction. [number] Default: 15
    -outSJfilterIntronMaxVsReadN	[s]	Maximum gap allowed for junctions supported by 1,2,3...N reads. [number] Default: 10.000.000
      
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
    'output'				=> undef,
    'nthreads'				=> 4,
    'fastqc'				=> "yes",
    'mapping'				=> "yes",
    'stranded'				=> "reversed",
    'fusionSearch'			=> "yes",
    'count'				=> "yes",
    'feature'				=> "gene_id",
    'uniq'				=> "no",
    'chimSegmentMin'			=> 15,
    'outSJfilterIntronMaxVsReadN'	=> 10000000,
    'chimJunctionOverhangMin'		=> 15,
    'species'				=> "HUMAN",
    'genome'				=> '/hpc/cog_bioinf/GENOMES',
    'fasta'				=> '/hpc/cog_bioinf/GENOMES',
    'fastqc_path'			=> '/hpc/cog_bioinf/common_scripts/FastQC/fastqc',
    'star_path'				=> '/hpc/local/CentOS6/cog_bioinf/STAR_2.3.0e/STAR',
    'samtools_path'			=> '/hpc/cog_bioinf/common_scripts/samtools-0.1.19/samtools',
    'sambamba_path'			=> '/hpc/local/CentOS6/cog_bioinf/bin/sambamba_v0.4.5',
    'qualimap_path'			=> '/hpc/cog_bioinf/common_scripts/qualimap_v0.7.1/qualimap',
    'picard_path'			=> '/hpc/cog_bioinf/common_scripts/picard-tools-1.98',
    'rnaseqc_path'			=> '/hpc/cog_bioinf/common_scripts/RNA-SeQC_v1.1.7.jar',
    'python_path'			=> '/hpc/local/CentOS6/cog_bioinf/Python-2.7.6/bin/python2.7',
    'gtf_file'				=> '/hpc/cog_bioinf/GENOMES',
    'merge_normalize_script'		=> '/hpc/cog_bioinf/common_scripts/RNA_seq_analysis/merge_normalize_count_tables.r',
    'refflat_file'			=> '/hpc/cog_bioinf/data/annelies/RNA_Seq'
);

die usage() if @ARGV == 0;
GetOptions (
    'h|help'				=> \$opt{help},
    'i|input=s'				=> \$opt{input},
    'o|output=s'			=> \$opt{output},
    'pe|nthreads=i'			=> \$opt{nthreads},
    'g|genome=s'			=> \$opt{genome},
    'fa|fasta=s'			=> \$opt{fasta},
    'stranded=s'			=> \$opt{stranded},
    'fqc|fastqc=s'			=> \$opt{fastqc},
    'm|mapping=s'			=> \$opt{mapping},
    'fusionSearch=s'			=> \$opt{fusionSearch},
    'c|count=s'				=> \$opt{count},
    'feature=s'				=> \$opt{feature},
    'uniq=s'				=> \$opt{uniq},
    'chimSegmentMin=i'			=> \$opt{chimSegmentMin},
    'outSJfilterIntronMaxVsReadN=s'	=> \$opt{outSJfilterIntronMaxVsReadN},
    'chimJunctionOverhangMin=i'		=> \$opt{chimJunctionOverhangMin},
    'sp|species=s'			=> \$opt{species}
) or die usage();

#check input paramaters
die usage() if $opt{help};
die "[ERROR] Number of threads must be at least 3!\n" if ($opt{nthreads} < 3);
die usage() unless $opt{input};
die usage() unless $opt{output};
die "[ERROR] Nothing to do!\n" if ( ($opt{fastqc} eq "no") && ($opt{mapping} eq "no") && ($opt{count} eq "no") );
die "[ERROR] Wrong option ($opt{fastqc}) given with -fastqc. Must be \"yes\" or \"no\".\n" if ( ($opt{fastqc} ne "no") && ($opt{fastqc} ne "yes") );
die "[ERROR] Wrong option ($opt{mapping}) given with -mapping. Must be \"yes\" or \"no\".\n" if ( ($opt{mapping} ne "no") && ($opt{mapping} ne "yes") );
die "[ERROR] Wrong option ($opt{stranded}) given with -stranded. Must be \"yes\", \"no\" or \"reversed\".\n" if ( ($opt{stranded} ne "no") && ($opt{stranded} ne "yes") && ($opt{stranded} ne "reversed") );
die "[ERROR] Wrong option ($opt{fusionSearch}) given with -fusionSearch. Must be \"yes\" or \"no\".\n" if ( ($opt{fusionSearch} ne "no") && ($opt{fusionSearch} ne "yes") );
die "[ERROR] Value given with -chimSegmentMin ($opt{chimSegmentMin}) must be higher than 0, when searching for fusion genes.\n" if ( ($opt{fusionSearch} eq "yes") && ($opt{chimSegmentMin} <= 0) );
die "[ERROR] Wrong option ($opt{count}) given with -count. Must be \"yes\" or \"no\".\n" if ( ($opt{count} ne "no") && ($opt{count} ne "yes") );
die "[ERROR] Wrong option ($opt{uniq}) given with -uniq. Must be \"yes\" or \"no\".\n" if ( ($opt{uniq} ne "no") && ($opt{uniq} ne "yes") );
die "[ERROR] Wrong option ($opt{feature}) given with -feature. Must be \"gene_id\", \"transcript_id\" or \"exon_id\".\n" if ( ($opt{feature} ne "gene_id") && ($opt{feature} ne "transcript_id") && ($opt{feature} ne "exon_id") );

my $SPECIES = uc $opt{species};
if ($SPECIES eq "HUMAN"){
    $opt{genome} .= '/Homo_sapiens.GRCh37.GATK.illumina/STAR';
    $opt{fasta} .= '/Homo_sapiens.GRCh37.GATK.illumina/STAR/hg19.fa';
    $opt{gtf_file} .= '/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.74.gtf';
    $opt{refflat_file} .= '/hg19.refFlat.gz';
} elsif ($SPECIES eq "RAT"){
    $opt{genome} .= '/rat_GATK_illumina_rnor_50/STAR';
    $opt{fasta} .= '/rat_GATK_illumina_rnor_50/STAR/rnor50.fa';
    $opt{gtf_file} .= '/rat_GATK_illumina_rnor_50/Rattus_norvegicus.Rnor_5.0.71.gtf';
    $opt{refflat_file} .= '/rnor50.refFlat.gz';
} elsif ($SPECIES eq "MOUSE"){
    $opt{genome} .= '/Mus_musculus_GRCm38_GATK_illumina_bwa075/STAR/mus_musculus_GRCm38';
    $opt{fasta} .= '/Mus_musculus_GRCm38_GATK_illumina_bwa075/STAR/mus_musculus_GRCm38/Mm_GRCm38_gatk_sorted.fa';
    $opt{gtf_file} .= '/Mus_musculus_GRCm38_GATK_illumina_bwa075/STAR/mus_musculus_GRCm38/Mus_musculus.GRCm38.70.gtf';
    $opt{refflat_file} .= '/Mus_musculus_GRCm38.refFlat.gz';
} elsif ($SPECIES eq "ZEBRAFISH"){
    $opt{genome} .= '/zfish9/STAR';
    $opt{fasta} .= '/zfish9/STAR/Zv9_66.fa';
    $opt{gtf_file} .= '/zfish9/Danio_rerio.Zv9.75.gtf';
    $opt{refflat_file} .= '/zfish9.refFlat.gz';
} else { 
    die "[ERROR] Wrong species ($SPECIES). Only HUMAN, RAT, MOUSE or ZEBRAFISH genomes are allowed.\n"
}



## ======================================
## Retrieve inputfiles
## ======================================

my @samplefiles;

foreach my $input ($opt{input}){
    open (FIND, "find -L $input -name '*_R1*.fastq.gz' |");

    while (my $f= <FIND>) {
	chomp $f;
	my $pattern = 'Undetermined_indices';
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
    my @path_parts = split("/", $f);
        
    my $base = basename($f);
    my @parts = split("_",$base);
    my $sample = join("_",@parts[0..$#parts-3]);
    my @input_parts = split("/","$opt{input}");
    my $run_name = $input_parts[$#input_parts];
    
    $rundir="$opt{output}/$run_name";
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
my $paired = 0;


foreach my $sample (keys %{$samples}) {
	
	print "Files for sample $sample:\n";
	my $job_id = "STAR_$sample\_".get_job_id();
	
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
	    $R2 =~ s/R1/R2/;
	    
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
	print STAR_SH "$opt{sambamba_path} sort -n -t $opt{nthreads} -o $rundir/$sample/mapping/$sample\_sorted_byName.bam $rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.out.rgAdded.bam\n";
	
	if ( $opt{uniq} eq 'yes' ){
	    print STAR_SH "$opt{sambamba_path} view -t $opt{nthreads} -f bam -F \"not secondary_alignment\" $rundir/$sample/mapping/$sample\_sorted.bam >$rundir/$sample/mapping/$sample\_sorted_uniq.bam\n";
	    print STAR_SH "java -Xmx2g  -jar $opt{picard_path}\/MarkDuplicates.jar INPUT=$rundir/$sample/mapping/$sample\_sorted_uniq.bam OUTPUT=$rundir/$sample/mapping/$sample\_sorted_uniq_markDup.bam METRICS_FILE=$rundir/$sample/mapping/$sample\_sorted_uniq_markDup_metrics.txt\n";
	    print STAR_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sorted_uniq.bam\n";
	    print STAR_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sorted_uniq_markDup.bam\n";
	    print STAR_SH "$opt{sambamba_path} flagstat $rundir/$sample/mapping/$sample\_sorted_uniq_markDup.bam\n";
	}
	#get stats from picard
	my $picard_strand = '';
	$picard_strand = 'FIRST_READ_TRANSCRIPTION_STRAND' if $opt{stranded} eq 'yes';
	$picard_strand = 'SECOND_READ_TRANSCRIPTION_STRAND' if $opt{stranded} eq 'reversed';
	$picard_strand = 'NONE' if $opt{stranded} eq 'no';
	print STAR_SH "java -jar $opt{picard_path}\/CollectRnaSeqMetrics.jar REF_FLAT=$opt{refflat_file} INPUT=$rundir/$sample/mapping/$sample\_sorted.bam OUTPUT=$rundir/$sample/mapping/$sample\_sorted_PicardRnaMetrics REFERENCE_SEQUENCE=$opt{fasta} STRAND_SPECIFICITY=$picard_strand\n";

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
	    
	    my $file = "$rundir/$sample/mapping/$sample\_sorted_byName.bam";
	    if ($opt{mapping} eq "no" && ! -e $file) {
		print HT "$opt{sambamba_path} sort -n -t $opt{nthreads} -o $rundir/$sample/mapping/$sample\_sorted_byName.bam $rundir/$sample/mapping/$sample\_sorted.bam\n";
	    }
	    
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
	    
	    print HT "$opt{sambamba_path} view $rundir/$sample/mapping/$sample\_sorted_byName.bam | $opt{python_path} -m HTSeq.scripts.count -m union -s $s_val -i $opt{feature} - $opt{gtf_file} > $rundir/$sample/read_counts/$sample\_htseq_counts.txt\n";

	    print HT "rm $rundir/$sample/mapping/$sample\_sorted_byName.bam*\n";
	    close HT;
	    if ( $opt{mapping} eq "no" ){
		print QSUB  "qsub -q veryshort -pe threaded $opt{nthreads} -o $rundir/$sample/logs -e $rundir/$sample/logs -N $HTSeqCount_job_id $rundir/$sample/jobs/$HTSeqCount_job_id.sh\n\n";
	    }else{
		print QSUB  "qsub -q veryshort -pe threaded $opt{nthreads} -o $rundir/$sample/logs -e $rundir/$sample/logs -N $HTSeqCount_job_id -hold_jid $job_id $rundir/$sample/jobs/$HTSeqCount_job_id.sh\n\n";
	    }
	}
}

if ( $opt{count} eq "yes"){
	#merge count tables of all samples and normalize the merged table
	my ($normTables_job_id) = ("normTables_".get_job_id());
	open MT,">$rundir/$normTables_job_id.sh";
	print MT "\#!/bin/sh\n\#\$ -S /bin/sh\n\n";
	print MT "uname -n > $rundir/$normTables_job_id.host\n";
	print MT "time R --save --args $rundir < $opt{merge_normalize_script}\n";
	close MT;
	my $hold_line = join ',',@hold_ids;
	print QSUB "qsub -q veryshort -o $rundir -e $rundir -N $normTables_job_id -hold_jid $hold_line $rundir/$normTables_job_id.sh\n\n";
}

close QSUB;
system "sh $mainJobID";


sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}

