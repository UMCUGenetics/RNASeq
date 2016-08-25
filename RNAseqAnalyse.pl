#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(tmpnam);
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use List::MoreUtils qw(uniq);
use Time::localtime;

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
    -sp|species                       [s]	Species of RNA-seq data [HUMAN/RAT/MOUSE/ZEBRAFISH/DOG/ARABIDOPSIS] Default: HUMAN
    -pe|nthreads                      [i]	Number of threads for processes run on execute nodes. [number] Default: 4
    -fqc|fastqc                       [s]	Perform FastQC? [yes/no] Default: yes
    -m|mapping                        [s]	Perform mapping with STAR? [yes/no] Default: yes
    -smallRNA                         [s]	Use STAR settings for small RNA dataset? [yes/no] Default: no
    -stranded                         [s]	Is the RNA-seq data from a strand-specific assay? [yes/no/reversed] Default: reversed
    -fusionSearch                     [s]	Want to detect fusion genes? [yes/no] Default: yes
    -c|count                          [s]	Want to count the mapped reads? [yes/no] Default: yes
    -id                               [s]	Which id attribute do you want to count? [gene_id/transcript_id/exon_id/all] Default: gene_id
    -n|normalize                      [s]	Want to normalize counted reads? [yes/no] Default: yes
    -rpkm                             [s]	Want to retrieve RPKMs? [yes/no] Default: yes
    -de                               [s]	Want to perform DE analysis? [yes/no] Default: no
    -design                           [s]	File with study design [/path/to/studydesign.txt] (see description below for required lay-out)
    -vc|variantcalling                [s]	Want to call variants? [yes/no] Default: no
    -bamqc                            [s]	Want to retrieve bam statistics? [yes/no] Default: yes
    -chimSegmentMin                   [i]	Minimum length of chimeric segment length. [number] Default: 15
    -chimJunctionOverhangMin          [i]	Minimum overhang for a chimeric junction. [number] Default: 15
    -outSJfilterIntronMaxVsReadN      [s]	Maximum gap allowed for junctions supported by 1,2,3...N reads. [number] Default: 10.000.000
    -debug                            [s]	Debug mode
    
    Lay-out design file:
     sample1	test
     sample2	control
     sample3	test
     sample4	test
     sample5	control
      
END
  exit;
}



## ======================================
## Get options
## ======================================

my %opt;
%opt = (
    'help'				=> undef,
    'debug'				=> undef,
    'input'				=> undef,
    'outputDir'				=> undef,
#     'mail'				=> undef,
    'nthreads'				=> 4,
    'fastqc'				=> "yes",
    'mapping'				=> "yes",
    'stranded'				=> "reversed",
    'fusionSearch'			=> "yes",
    'count'				=> "yes",
    'merge'				=> "yes",
    'normalize'				=> "yes",
    'rpkm'				=> "yes",
    'bamqc'				=> "yes",
    'de'				=> "no",
    'design'				=> undef,
    'variantcalling'			=> "no",
    'id'				=> "gene_id",
    'uniq'				=> "no",
    'chimSegmentMin'			=> 15,
    'outSJfilterIntronMaxVsReadN'	=> 10000000,
    'chimJunctionOverhangMin'		=> 15,
    'smallRNA'				=> "no",
    'species'				=> "HUMAN",
    'genome'				=> '/hpc/cog_bioinf/GENOMES/STAR',
    'fasta'				=> undef,
    'intervallist'			=> undef,
    'gtf_file'				=> undef,
    'refflat_file'			=> undef,
    'genesizes_file'			=> undef,
    'fastqc_path'			=> '/hpc/local/CentOS7/cog_bioinf/FastQ-v0.11.4/fastqc',
    'star_path'				=> '/hpc/local/CentOS7/cog_bioinf/STAR-STAR_2.4.2a/source/STAR',
    'sambamba_path'			=> '/hpc/local/CentOS7/cog_bioinf/sambamba_v0.5.8/sambamba',
    'picard_path'			=> '/hpc/local/CentOS7/cog_bioinf/picard-tools-1.141/picard.jar',
    'gatk_path'				=> '/hpc/local/CentOS7/cog_bioinf/GenomeAnalysisTK-3.4-46',
    'bamstats_path'			=> '/hpc/local/CentOS7/cog_bioinf/bamMetrics/bamMetrics.pl',
    'snpEff_path'			=> '/hpc/local/CentOS7/cog_bioinf/snpEff_v4.1h',
    'igvtools_path'			=> '/hpc/local/CentOS7/cog_bioinf/igvtools-2.3.60/igvtools',
    'dbnsfp_path'			=> '/hpc/cog_bioinf/common_dbs/dbNSFP/dbNSFPv2.9/dbNSFP2.9.txt.gz',
);

die usage() if @ARGV == 0;
GetOptions (
    'h|help'				=> \$opt{help},
    'debug'				=> \$opt{debug},
    'i|input=s@'			=> \$opt{input},
    'o|outputDir=s'			=> \$opt{outputDir},
#     'mail=s'				=> \$opt{mail},
    'pe|nthreads=i'			=> \$opt{nthreads},
    'g|genome=s'			=> \$opt{genome},
    'fa|fasta=s'			=> \$opt{fasta},
    'stranded=s'			=> \$opt{stranded},
    'fqc|fastqc=s'			=> \$opt{fastqc},
    'm|mapping=s'			=> \$opt{mapping},
    'fusionSearch=s'			=> \$opt{fusionSearch},
    'c|count=s'				=> \$opt{count},
    'merge=s'				=> \$opt{merge},
    'n|normalize=s'			=> \$opt{normalize},
    'rpkm=s'				=> \$opt{rpkm},
    'de=s'				=> \$opt{de},
    'design=s'				=> \$opt{design},
    'vc|variantcalling=s'		=> \$opt{variantcalling},
    'smallRNA=s'			=> \$opt{smallRNA},
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
die "[ERROR] Nothing to do!\n" if ( ($opt{fastqc} eq "no") && ($opt{mapping} eq "no") && ($opt{count} eq "no") && ($opt{normalize} eq "no") && ($opt{rpkm} eq "no") && ($opt{bamqc} eq "no") && ($opt{de} eq "no") && ($opt{vc} eq "no") );
die "[ERROR] Wrong option ($opt{fastqc}) given with -fastqc. Must be \"yes\" or \"no\".\n" if ( ($opt{fastqc} ne "no") && ($opt{fastqc} ne "yes") );
die "[ERROR] Wrong option ($opt{mapping}) given with -mapping. Must be \"yes\" or \"no\".\n" if ( ($opt{mapping} ne "no") && ($opt{mapping} ne "yes") );
die "[ERROR] Wrong option ($opt{stranded}) given with -stranded. Must be \"yes\", \"no\" or \"reversed\".\n" if ( ($opt{stranded} ne "no") && ($opt{stranded} ne "yes") && ($opt{stranded} ne "reversed") );
die "[ERROR] Wrong option ($opt{fusionSearch}) given with -fusionSearch. Must be \"yes\" or \"no\".\n" if ( ($opt{fusionSearch} ne "no") && ($opt{fusionSearch} ne "yes") );
die "[ERROR] Value given with -chimSegmentMin ($opt{chimSegmentMin}) must be higher than 0, when searching for fusion genes.\n" if ( ($opt{fusionSearch} eq "yes") && ($opt{chimSegmentMin} <= 0) );
die "[ERROR] Wrong option ($opt{count}) given with -count. Must be \"yes\" or \"no\".\n" if ( ($opt{count} ne "no") && ($opt{count} ne "yes") );
die "[ERROR] Wrong option ($opt{normalize}) given with -normalize. Must be \"yes\" or \"no\".\n" if ( ($opt{normalize} ne "no") && ($opt{normalize} ne "yes") );
die "[ERROR] Wrong option ($opt{rpkm}) given with -rpkm. Must be \"yes\" or \"no\".\n" if ( ($opt{rpkm} ne "no") && ($opt{rpkm} ne "yes") );
die "[ERROR] Wrong option ($opt{uniq}) given with -uniq. Must be \"yes\" or \"no\".\n" if ( ($opt{uniq} ne "no") && ($opt{uniq} ne "yes") );
die "[ERROR] Wrong option ($opt{id}) given with -id. Must be \"gene_id\", \"transcript_id\", \"exon_id\" or \"all\".\n" if ( ($opt{id} ne "gene_id") && ($opt{id} ne "transcript_id") && ($opt{id} ne "exon_id") && ($opt{id} ne "all"));
die "[ERROR] Wrong option ($opt{variantcalling}) given with -variantcalling. Must be \"yes\" or \"no\".\n" if ( ($opt{variantcalling} ne "no") && ($opt{variantcalling} ne "yes") );
die "[ERROR] Wrong option ($opt{smallRNA}) given with -variantcalling. Must be \"yes\" or \"no\".\n" if ( ($opt{smallRNA} ne "no") && ($opt{smallRNA} ne "yes") );
die "[ERROR] Must give design (-design) with -de option.\n" if ( ($opt{de} eq "yes") && (! $opt{design}) );
die "[ERROR] Design file doesn't exist!\n" if ( ($opt{de} eq "yes") && (! -e $opt{design}) );
die "[ERROR] Design file is not in right format. See -h for lay-out.\n" if ( ($opt{de} eq "yes") && (checkDesign($opt{design}) eq "wrong") );

my $SPECIES = uc $opt{species};
my %annotationDB = ( 'HUMAN' => 'org.Hs.eg.db', 'RAT' => 'org.Rn.eg.db', 'MOUSE' => 'org.Mm.eg.db', 'ZEBRAFISH' => 'org.Dr.eg.db', 'DOG' => 'org.Cf.eg.db', 'ARABIDOPSIS' => 'org.At.tair.db' );
my @knownSites = qw(/hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/1000G_phase1.indels.b37.vcf /hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/dbsnp_137.b37.vcf /hpc/cog_bioinf/common_scripts/GATK_v2.7/bundle/Mills_and_1000G_gold_standard.indels.b37.vcf);

if ($SPECIES eq "HUMAN"){
    $opt{genome} .= '/Homo_sapiens.GRCh37';
    $opt{fasta} = $opt{genome}.'/Homo_sapiens.GRCh37.GATK.illumina.fa';
    $opt{gtf_file} = $opt{genome}.'/Homo_sapiens.GRCh37.74.gtf';
    $opt{refflat_file} = $opt{genome}.'/hg19.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/Homo_sapiens.GRCh37.GATK.illumina.rRNA.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Homo_sapiens.GRCh37.74_exon_gene_sizes.txt';
    $opt{annotate_db} = 'GRCh37.74';
} elsif ($SPECIES eq "RAT"){
    $opt{genome} .= '/Rattus_norvegicus.Rnor50';
    $opt{fasta} = $opt{genome}.'/Rn_Rn05_ill_gatk_sorted.fa';
    $opt{gtf_file} = $opt{genome}.'/Rattus_norvegicus.Rnor_5.0.71.gtf';
    $opt{refflat_file} = $opt{genome}.'/rnor50.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/rno5_rRNA_intervallist.txt';
    $opt{genesizes_file} = $opt{genome}.'/Rattus_norvegicus.Rnor_5.0.71_exon_gene_sizes.txt';
    $opt{annotate_db} = 'Rnor_5.0.71';
} elsif ($SPECIES eq "MOUSE"){
    $opt{genome} .= '/Mus_musculus.GRCm38';
    $opt{fasta} = $opt{genome}.'/Mm_GRCm38_gatk_sorted.fa';
    $opt{gtf_file} = $opt{genome}.'/Mus_musculus.GRCm38.70.gtf';
    $opt{refflat_file} = $opt{genome}.'/Mus_musculus_GRCm38.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/Mus_musculus_GRCm38.rRNA.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Mus_musculus.GRCm38.70_exon_gene_sizes.txt';
    $opt{annotate_db} = 'GRCm38.71';
} elsif ($SPECIES eq "ZEBRAFISH"){
    $opt{genome} .= '/Danio_rerio.Zv9';
    $opt{fasta} = $opt{genome}.'/Zv9_66.fa';
    $opt{gtf_file} = $opt{genome}.'/Danio_rerio.Zv9.75.gtf';
    $opt{refflat_file} = $opt{genome}.'/zfish9.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/Danio_rerio.Zv9.75.rRNA.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Danio_rerio.Zv9.75_exon_gene_sizes.txt';
    $opt{annotate_db} = 'Zv9.75';
} elsif ($SPECIES eq "DOG"){
    $opt{genome} .= '/Canis_familiaris.CanFam31';
    $opt{fasta} = $opt{genome}.'/Cf3_ens71_GATK.fa';
    $opt{gtf_file} = $opt{genome}.'/Canis_familiaris.CanFam3.1.75.gtf';
    $opt{refflat_file} = $opt{genome}.'/CanFam3.1.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/CanFam3.1_rRNA_genes.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Canis_familiaris.CanFam3.1.75_exon_gene_sizes.txt';
    $opt{annotate_db} = 'CanFam3.1.71';
} elsif ($SPECIES eq "ARABIDOPSIS"){
    $opt{genome} .= '/Arabidopsis_thaliana.TAIR10';
    $opt{fasta} = $opt{genome}.'/Arabidopsis_thaliana.TAIR10.30.fa';
    $opt{gtf_file} = $opt{genome}.'/Arabidopsis_thaliana.TAIR10.30.gtf';
    $opt{refflat_file} = $opt{genome}.'/Arabidopsis_thaliana.TAIR10.30.refFlat.gz';
    $opt{intervallist} = $opt{genome}.'/Arabidopsis_thaliana.TAIR10.30.rRNA.intervallist';
    $opt{genesizes_file} = $opt{genome}.'/Arabidopsis_thaliana.TAIR10.30_exon_gene_sizes.txt';
    $opt{annotate_db} = 'athalianaTair10';
}else {
    die "[ERROR] Wrong species ($SPECIES). Only HUMAN, RAT, MOUSE, DOG, ZEBRAFISH or ARABIDOPSIS genomes are allowed.\n";
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
    my @fastqfiles = glob($fastqdir."/*{/,}*_R1_*.fastq.gz");
    foreach my $fastq (@fastqfiles){
	my $pattern = 'Undertermined';
	push @samplefiles, $fastq unless $fastq =~ /$pattern/;
    }
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
my @final_bams;

foreach my $f (@samplefiles){
        
    my $base = basename($f);
    my @parts = split("_",$base);
    my $sample = $parts[0];
    my @input_parts = split("/","$opt{input}");
    
    $rundir=$opt{outputDir};
    push(@{$samples->{$sample}}, $f);
        
    if(! -e "$rundir"){
	make_path($rundir) or die "Couldn't create run directory: $rundir\n";
    }
    if(! -e "$rundir/read_counts" && $opt{count} ne 'no'){
	mkdir("$rundir/read_counts") or die "Couldn't create readcount directory: $rundir/read_counts\n";
    }
    if(! -e "$rundir/logs"){
	mkdir("$rundir/logs") or die "Couldn't create logs directory: $rundir/logs\n";
    }
    if(! -e "$rundir/jobs"){
	mkdir("$rundir/jobs") or die "Couldn't create jobs directory: $rundir/jobs\n";
    }
    if(! -e "$rundir/$sample"){
	mkdir("$rundir/$sample") or die "Couldn't create sample directory: $rundir/$sample\n";
    }
    if(! -e "$rundir/$sample/fastqc" && $opt{fastqc} ne 'no'){
        mkdir("$rundir/$sample/fastqc") or die "Couldn't create fastqc directory: $rundir/$sample/fastqc\n";
    }
    if(! -e "$rundir/$sample/mapping" && $opt{mapping} ne 'no'){
        mkdir("$rundir/$sample/mapping") or die "Couldn't create mapping directory: $rundir/$sample/mapping\n";
    }
    if(! -e "$rundir/$sample/read_counts" && $opt{count} ne 'no'){
	mkdir("$rundir/$sample/read_counts") or die "Couldn't create readcounts directory: $rundir/$sample/read_counts\n";
    }
#     if(! -e "$rundir/variantCalling" && $opt{variantcalling} ne 'no'){
# 	mkdir("$rundir/variantCalling") or die "Couldn't create variants directory: $rundir/variantCalling\n";
#     }
    if(! -e "$rundir/$sample/jobs"){
        mkdir("$rundir/$sample/jobs") or die "Couldn't create jobs directory: $rundir/$sample/jobs\n";
    }
    if(! -e "$rundir/$sample/logs"){
        mkdir("$rundir/$sample/logs") or die "Couldn't create logs directory: $rundir/$sample/logs\n";
    }
    if(! -e "$rundir/DEanalysis" && $opt{de} eq 'yes'){
	mkdir("$rundir/DEanalysis") or die "Couldn't create DE analysis directory: $rundir/DEanalysis\n";
	system("cp $opt{design} $rundir/DEanalysis/");
    }
}



## ======================================
## Create settings file
## ======================================

#check pipeline version
my $PIPELINE_PATH = dirname(abs_path($0));
chdir $PIPELINE_PATH;
my $PIPELINE_VERSION = qx(git describe --always --tag | cut -d '-' -f 1);

#copy readme to run directory
system("cp $PIPELINE_PATH/README.md $rundir");

#write settings info to settings.txt
my $timestamp = timestamp();
open SETTINGS, ">$rundir/settings_$timestamp.txt" or die "Cannot open settings file!\n";
my $datestring = localtime();
print SETTINGS "$datestring\n\n";
print SETTINGS "VERSION=$PIPELINE_VERSION\n";
print SETTINGS "SPECIES=$opt{species}\n";
print SETTINGS "RUNNAME=".basename($rundir)."\n";
print SETTINGS "FASTQC=$opt{fastqc}\n";
print SETTINGS "MAPPING=$opt{mapping}\n";
print SETTINGS "FUSION=$opt{fusionSearch}\n";
print SETTINGS "COUNTING=$opt{count}\n";
print SETTINGS "RPKM=$opt{rpkm}\n";
print SETTINGS "DE=$opt{de}\n";
print SETTINGS "BAMQC=$opt{bamqc}\n";
print SETTINGS "VARIANTCALLING=$opt{variantcalling}\n";



## ======================================
## Create jobs
## ======================================

my $mainJobID = "$rundir/jobs/".get_job_id()."_qsub.sh";
open QSUB, ">$mainJobID";
print QSUB "\#!/bin/bash\n\#\$ \-o $rundir/logs\n\#\$ \-e $rundir/logs\n\n";

#cleanup script
my ($clean_job_id) = ("clean_".get_job_id());
open CL, ">$rundir/jobs/$clean_job_id.sh";
print CL "\#!/bin/bash\n\n";
print CL "uname -n > $rundir/logs/$clean_job_id.host\n";

my @hold_ids=();
my @hold_mapping_ids=();
my $paired = 0;
my $runname = basename( $rundir );

foreach my $sample (keys %{$samples}) {
	print "\nFiles for sample $sample:\n";
	print SETTINGS "\nFiles for sample $sample:\n";
	print QSUB "\n\n###Submissions for sample $sample:\n";
	my $job_id = "MAPPING_$sample\_".get_job_id();
	push @hold_mapping_ids, $job_id;
	#create bash script for STAR submission of this sample
	open MAPPING_SH,">$rundir/$sample/jobs/$job_id.sh" or die "Couldn't create $rundir/$sample/jobs/$job_id.sh\n";
	print MAPPING_SH "\#!/bin/bash\n\n";
	print MAPPING_SH "uname -n > $rundir/$sample/logs/$sample.host\n";
	print MAPPING_SH "echo \"mapping pair\t\" `date` >> $rundir/$sample/logs/$sample.host\n";
	print MAPPING_SH "mkdir -p $rundir/$sample/mapping/tmp/$job_id/\n";
	print MAPPING_SH "cd $rundir/$sample/mapping/tmp/$job_id/\n\n";
	#Create STAR command
	my $star_command = "$opt{star_path} --genomeDir $opt{genome} --runThreadN $opt{nthreads} --outFileNamePrefix $sample\_ --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat";
	if ($opt{stranded} eq "no"){
	    $star_command .= " --outSAMstrandField intronMotif";
	}
	if (defined $opt{outSJfilterIntronMaxVsReadN}) {
	    $star_command .= " --outSJfilterIntronMaxVsReadN " . $opt{outSJfilterIntronMaxVsReadN};
	}
	if (defined $opt{chimJunctionOverhangMin}) {
	    $star_command .= " --chimJunctionOverhangMin " . $opt{chimJunctionOverhangMin};
	}
	if ($opt{smallRNA} eq "yes"){
	    $star_command .= " --outFilterMultimapNmax 50 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.3 --outFilterMatchNmin 16 --outFilterMatchNminOverLread 0 --outFilterScoreMinOverLread 0 --clip3pAdapterSeq TGGAATTCTCGGGTGCCAAGG --clip3pAdapterMMp 0.1 --clip5pNbases 0 --alignIntronMax 1.0 --alignEndsType Local";
	}
	if ($opt{fusionSearch} eq "yes"){
	    $star_command .= " --chimSegmentMin $opt{chimSegmentMin}";
	}
	elsif ($opt{fusionSearch} eq "no"){
	    $star_command .= " --chimSegmentMin 0";
	}
	if ($opt{variantcalling} eq "yes"){
	    $star_command .= " --twopassMode Basic";
	}
	
	my ($ID,$PU,$LB,$SM,$PL)='';
	
	#Find Fastq files
	my $R1_fastqs = "";
	my $R2_fastqs = "";
	my $rest_fastqs = "";
	foreach my $fastq (@{$samples->{$sample}}){
	    my $R1 = $fastq;
	    my $R2 = $fastq;
	    $R2 =~ s/_R1_/_R2_/;
	    
	    my $fastqname = basename($R1);
	    my @cols = split("_",$fastqname);
	    $ID = join("_",@cols[0..2]);
	    $PU = $cols[1];
	    $LB = $cols[2];
	    $SM = $cols[0];
	    $PL = 'ILLUMINA';
	    
	    #determine paired or single end
	    if (-e $R2) {
		$paired = 1;
		$R1_fastqs .= "$R1,";
		$R2_fastqs .= "$R2,";
	    
		print "Pair found:\n";
		print $R1," = R1\n";
		print $R2," = R2\n\n";
		print SETTINGS "Pair found:\n";
		print SETTINGS $R1," = R1\n";
		print SETTINGS $R2," = R2\n\n";
	    }else{
		$rest_fastqs .= "$R1,";
		print "Single tag found\n";
		print $R1," = R1\n";
		print SETTINGS "Single tag found\n";
		print SETTINGS $R1," = R1\n";
	    }
	    
	    #FASTQC
	    if ($opt{fastqc} eq "yes") {
		my ($FastQC1_job_id, $FastQC2_job_id) = ("FastQC1_".get_job_id(), "FastQC2_".get_job_id());
		
		
		open FASTQC,">$rundir/$sample/jobs/$FastQC1_job_id.sh";
		print FASTQC "\#!/bin/bash\n\n";
		print FASTQC "cd $rundir/$sample/fastqc\n\n";
		print FASTQC "uname -n > ../logs/$FastQC1_job_id.host\n";
		print FASTQC "echo \"FastQC\t\" `date` >> ../logs/FastQC_$sample.host\n";
		print FASTQC "$opt{fastqc_path} $R1 -o $rundir/$sample/fastqc\n";
		print QSUB "##FastQC\nqsub -l h_rt=02:00:00 -o $rundir/$sample/logs -e $rundir/$sample/logs -R yes -N $FastQC1_job_id $rundir/$sample/jobs/$FastQC1_job_id.sh\n";
		close FASTQC;
		
		
		if ($paired == 1) {
		    open FASTQC2,">$rundir/$sample/jobs/$FastQC2_job_id.sh";
		    print FASTQC2 "\#!/bin/bash\n\n";
		    print FASTQC2 "cd $rundir/$sample/fastqc\n\n";
		    print FASTQC2 "uname -n > ../logs/$FastQC2_job_id.host\n";
		    print FASTQC2 "echo \"FastQC\t\" `date` >> ../logs/FastQC2_$sample.host\n";
		    print FASTQC2 "$opt{fastqc_path} $R2 -o $rundir/$sample/fastqc\n";
		    print QSUB "##FastQC\nqsub -l h_rt=02:00:00 -o $rundir/$sample/logs -e $rundir/$sample/logs -R yes -N $FastQC2_job_id $rundir/$sample/jobs/$FastQC2_job_id.sh\n";
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
	print MAPPING_SH $star_command;
	
	#add read groups
	print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{picard_path} AddOrReplaceReadGroups INPUT=$rundir/$sample/mapping/tmp/$job_id/$sample\_Aligned.sortedByCoord.out.bam OUTPUT=$rundir/$sample/mapping/$sample\_sort.bam RGID=$ID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM\n";
	print MAPPING_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort.bam\n";
	my $finalbam = "$rundir/$sample/mapping/$sample\_sort.bam";
	
	if ( $opt{variantcalling} eq 'yes'){
	    #mark duplicates
	    print MAPPING_SH "\n###Mark Duplicates\n";
	    print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{picard_path} MarkDuplicates INPUT=$rundir/$sample/mapping/$sample\_sort.bam OUTPUT=$rundir/$sample/mapping/$sample\_sort_dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$rundir/$sample/mapping/$sample\_markDup_metrics.txt\n";
	    print MAPPING_SH "$opt{sambamba_path} flagstat -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort_dedup.bam > $rundir/$sample/mapping/$sample\_sort_dedup.flagstat\n";
	    print MAPPING_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort_dedup.bam\n";
	    #splitNCigarReads
	    print MAPPING_SH "\n###SplitNTrim\n";
	    print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T SplitNCigarReads -R $opt{fasta} -I $rundir/$sample/mapping/$sample\_sort_dedup.bam -o $rundir/$sample/mapping/$sample\_sort_dedup_splitN.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS\n";
	    print MAPPING_SH "$opt{sambamba_path} flagstat -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort_dedup_splitN.bam > $rundir/$sample/mapping/$sample\_sort_dedup_splitN.flagstat\n";
	    print MAPPING_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort_dedup_splitN.bam\n";
	    #realigner target creator
	    print MAPPING_SH "\n###Indel Realignment\n";
	    print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $opt{fasta} -I $rundir/$sample/mapping/$sample\_sort_dedup_splitN.bam -o $rundir/$sample/mapping/$sample\_target.intervals\n";
	    #indel realignment
	    print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T IndelRealigner -R $opt{fasta} -I $rundir/$sample/mapping/$sample\_sort_dedup_splitN.bam -targetIntervals $rundir/$sample/mapping/$sample\_target.intervals -o $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned.bam\n";
	    print MAPPING_SH "$opt{sambamba_path} flagstat -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned.bam > $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned.flagstat\n";
	    print MAPPING_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned.bam\n";
	    $finalbam = "$rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned.bam";
	    #baserecalibrator
	    if ( $SPECIES eq 'HUMAN' ){
		print MAPPING_SH "\n###Base Recalibration\n";
		print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T BaseRecalibrator -R $opt{fasta} -I $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned.bam -knownSites " . join(" -knownSites ", @knownSites) . " -o $rundir/$sample/mapping/$sample\_recal_data.table\n";
		print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T BaseRecalibrator -R $opt{fasta} -I $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned.bam -knownSites " . join(" -knownSites ", @knownSites) . " -BQSR $rundir/$sample/mapping/$sample\_recal_data.table -o $rundir/$sample/mapping/$sample\_post_recal_data.table\n";
		print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $opt{fasta} -before $rundir/$sample/mapping/$sample\_recal_data.table -after $rundir/$sample/mapping/$sample\_post_recal_data.table -plots $rundir/$sample/mapping/$sample\_recalibration_plots.pdf\n";
		print MAPPING_SH "java -Xmx32G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T PrintReads -R $opt{fasta} -I $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned.bam -BQSR $rundir/$sample/mapping/$sample\_recal_data.table -o $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned_recal.bam\n";
		print MAPPING_SH "$opt{sambamba_path} flagstat -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned_recal.bam > $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned_recal.flagstat\n";
		print MAPPING_SH "$opt{sambamba_path} index -t $opt{nthreads} $rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned_recal.bam\n";
		$finalbam = "$rundir/$sample/mapping/$sample\_sort_dedup_splitN_realigned_recal.bam";
	    }
	}
	push @final_bams, $finalbam;
	
	#index final bam
	#print MAPPING_SH "$opt{sambamba_path} index -t $opt{nthreads} $finalbam\n";

	#remove and move files
	print CL "rm $rundir/$sample/mapping/$sample\_*Aligned.sortedByCoord.out.bam\n";
	print MAPPING_SH "mv $rundir/$sample/mapping/tmp/$job_id/* $rundir/$sample/mapping/\n";
	print CL "rm -r $rundir/$sample/mapping/tmp\n";
	print CL "rm -r $rundir/$sample/mapping/$sample\__STARtmp\n";
	
	close MAPPING_SH;
	
	print QSUB "\n";
	
	if ($opt{mapping} eq "yes") {
	    print QSUB "##MAPPING\nqsub -l h_rt=08:00:00,h_vmem=100G -pe threaded $opt{nthreads} -o $rundir/$sample/logs -e $rundir/$sample/logs -R yes -N $job_id $rundir/$sample/jobs/$job_id.sh\n\n";
	}
	
	#get read counts of mapped reads
	if ($opt{count} eq "yes") {
	    my ($HTSeqCount_job_id) = ("htseqCount_".get_job_id());
	    push (@hold_ids, $HTSeqCount_job_id);
	    open HT,">$rundir/$sample/jobs/$HTSeqCount_job_id.sh";
	    print HT "\#!/bin/bash\n\n";
	    print HT "uname -n > $rundir/$sample/logs/htseqCount_$sample.host\n";
	    print HT "module load python/2.7.10\n";
	    
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
	    
	    if ( $opt{id} eq "all" ){
		print HT "$opt{sambamba_path} view $rundir/$sample/mapping/$sample\_sort.bam | python -m HTSeq.scripts.count -m union -r pos -s $s_val -i gene_id - $opt{gtf_file} > $rundir/$sample/read_counts/$sample\_htseq_gene_counts.txt\n";
		print HT "$opt{sambamba_path} view $rundir/$sample/mapping/$sample\_sort.bam | python -m HTSeq.scripts.count -m union -r pos -s $s_val -i exon_id - $opt{gtf_file} > $rundir/$sample/read_counts/$sample\_htseq_exon_counts.txt\n";
		print HT "$opt{sambamba_path} view $rundir/$sample/mapping/$sample\_sort.bam | python -m HTSeq.scripts.count -m union -r pos -s $s_val -i transcript_id - $opt{gtf_file} > $rundir/$sample/read_counts/$sample\_htseq_transcript_counts.txt\n";
		print HT "module unload python/2.7.10\n";
		close HT;
		if ( $opt{mapping} eq "no" ){
		    print QSUB  "##HTSeq-count\nqsub -l h_rt=08:00:00 -o $rundir/$sample/logs -e $rundir/$sample/logs -R yes -N $HTSeqCount_job_id $rundir/$sample/jobs/$HTSeqCount_job_id.sh\n\n";
		}else{
		    print QSUB  "##HTSeq-count\nqsub -l h_rt=08:00:00 -o $rundir/$sample/logs -e $rundir/$sample/logs -R yes -N $HTSeqCount_job_id -hold_jid $job_id $rundir/$sample/jobs/$HTSeqCount_job_id.sh\n\n";
		}
	}
	    else{
		print HT "$opt{sambamba_path} view $rundir/$sample/mapping/$sample\_sort.bam | python -m HTSeq.scripts.count -m union -r pos -s $s_val -i $opt{id} - $opt{gtf_file} > $rundir/$sample/read_counts/$sample\_htseq_counts.txt\n";
		print HT "module unload python/2.7.10\n";
		close HT;
		if ( $opt{mapping} eq "no" ){
		    print QSUB  "##HTSeq-count\nqsub -l h_rt=08:00:00 -o $rundir/$sample/logs -e $rundir/$sample/logs -R yes -N $HTSeqCount_job_id $rundir/$sample/jobs/$HTSeqCount_job_id.sh\n\n";
		}else{
		    print QSUB  "##HTSeq-count\nqsub -l h_rt=08:00:00 -o $rundir/$sample/logs -e $rundir/$sample/logs -R yes -N $HTSeqCount_job_id -hold_jid $job_id $rundir/$sample/jobs/$HTSeqCount_job_id.sh\n\n";
		}
	    }
	}
}

if ( $opt{variantcalling} eq "yes" ){
    my ($variantcalling_job_id) = ("variantCalling_".get_job_id());
    open VC, ">$rundir/jobs/$variantcalling_job_id.sh";
    print VC "\#!/bin/bash\n\n";
    print VC "uname -n > $rundir/logs/$variantcalling_job_id.host\n";
    
    #variant calling
    my $rawvariants_vcf = "$rundir/".basename($rundir).".raw_variants.vcf";
    print VC "if [ ! -s $rawvariants_vcf ]\nthen\n\tjava -Xmx10G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T HaplotypeCaller -R $opt{fasta} -I ". join(" -I ",@final_bams) ." -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $rawvariants_vcf\nfi\n\n";
    
    #variant filtering
    my $filtered_vcf = "$rundir/".basename($rundir).".filtered_variants.vcf";
    print VC "if [ -s $rawvariants_vcf ] && [ ! -s $filtered_vcf ]\nthen\n\tjava -Xmx10G -Djava.io.tmpdir=\$TMPDIR -jar $opt{gatk_path}\/GenomeAnalysisTK.jar -T VariantFiltration -R $opt{fasta} -V $rawvariants_vcf -o $filtered_vcf -window 10 -cluster 3 -filterName LowQualityDepth --filterExpression \"QD < 2.0\" -filterName MappingQuality --filterExpression \"MQ < 40.0\" -filterName StrandBias --filterExpression \"FS > 60.0\" -filterName HaplotypeScoreHigh --filterExpression \"HaplotypeScore > 13.0\" -filterName MQRankSumLow --filterExpression \"MQRankSum < -12.5\" -filterName ReadPosRankSumLow --filterExpression \"ReadPosRankSum < -8.0\" -filterName HardToValidate --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" -filterName LowCoverage --filterExpression \"DP < 5\" -filterName VeryLowQual --filterExpression \"QUAL < 30\" -filterName LowQual --filterExpression \"QUAL >= 30.0 && QUAL < 50.0\" -filterName SOR --filterExpression \"SOR > 4.0\"\n";
    print VC "else\n\techo \"ERROR: $rawvariants_vcf does not exist.\"\nfi\n\n";
    
    #variant annotation
    #snpEff
    my $snpeff_vcf = "$rundir/".basename($rundir).".filtered_variants_snpEff.vcf";
    print VC "if [ -s $filtered_vcf ] && [ ! -s $snpeff_vcf ]\nthen\n\tjava -Xmx15G -Djava.io.tmpdir=\$TMPDIR -jar ".$opt{snpEff_path}."/snpEff.jar -c ".$opt{snpEff_path}."/snpEff.config $opt{annotate_db} -v $filtered_vcf -hgvs -lof -no-downstream -no-upstream -no-intergenic > $snpeff_vcf\n";
    #igvtools index
    print VC "\t$opt{igvtools_path} index $snpeff_vcf\n\trm igv.log\n\n";
    print VC "else\n\techo \"ERROR: $filtered_vcf does not exist.\"\nfi\n\n";
    
    
    if ( $SPECIES eq 'HUMAN' ){
	#SnpSift
	my $snpsift_vcf = "$rundir/".basename($rundir).".filtered_variants_snpEff_snpSift.vcf";
	my $snpsift_fields = "hg38_chr,hg38_pos,genename,Uniprot_acc,Uniprot_id,Uniprot_aapos,Interpro_domain,cds_strand,refcodon,SLR_test_statistic,codonpos,fold-degenerate,Ancestral_allele,Ensembl_geneid,Ensembl_transcriptid,aapos,aapos_SIFT,aapos_FATHMM,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationAssessor_score,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_score,FATHMM_rankscore,FATHMM_pred,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,Reliability_index,VEST3_score,VEST3_rankscore,PROVEAN_score,PROVEAN_converted_rankscore,PROVEAN_pred,CADD_raw,CADD_raw_rankscore,CADD_phred,GERP++_NR,GERP++_RS,GERP++_RS_rankscore,phyloP46way_primate,phyloP46way_primate_rankscore,phyloP46way_placental,phyloP46way_placental_rankscore,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore,phastCons46way_primate,phastCons46way_primate_rankscore,phastCons46way_placental,phastCons46way_placental_rankscore,phastCons100way_vertebrate,phastCons100way_vertebrate_rankscore,SiPhy_29way_pi,SiPhy_29way_logOdds,SiPhy_29way_logOdds_rankscore,LRT_Omega,UniSNP_ids,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF,ARIC5606_AA_AC,ARIC5606_AA_AF,ARIC5606_EA_AC,ARIC5606_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,clinvar_rs,clinvar_clnsig,clinvar_trait,COSMIC_ID,COSMIC_CNT";
	print VC "if [ -s $snpeff_vcf ] && [ ! -s $snpsift_vcf ]\nthen\n\tjava -Xmx15G -Djava.io.tmpdir=\$TMPDIR -jar ".$opt{snpEff_path}."/SnpSift.jar dbnsfp -v -f $snpsift_fields -db $opt{dbnsfp_path} $snpeff_vcf > $snpsift_vcf\n";
	#igvtools index
	print VC "\t$opt{igvtools_path} index $snpsift_vcf\n\trm igv.log\n\n";
	print VC "else\n\techo \"ERROR: $snpeff_vcf does not exist.\"\nfi\n\n";
	print VC "if [ -s $snpsift_vcf ]\nthen\n\trm $snpeff_vcf $snpeff_vcf.idx\nfi\n\n";
    }
    my $hold_line = join ',',@hold_mapping_ids;
    print QSUB "\n##Variant Calling\nqsub -l h_vmem=20G -pe threaded 4 -l h_rt=168:0:0 -o $rundir/logs -e $rundir/logs -R yes -N $variantcalling_job_id -hold_jid $hold_line $rundir/jobs/$variantcalling_job_id.sh\n\n";
}

if ( $opt{count} eq "yes" ){
    #merge count tables of all samples
    my ($mergeTables_job_id) = ("mergeTables_".get_job_id());
    open MT,">$rundir/jobs/$mergeTables_job_id.sh";
    print MT "\#!/bin/bash\n\n";
    print MT "uname -n > $rundir/logs/$mergeTables_job_id.host\n";
    #print MT "export MODULEPATH=/hpc/local/CentOS7/cog_bioinf/modules:\${MODULEPATH}\n";
    print MT "module load R/3.2.2\n";
    mergeRscript($rundir);
    print MT "time R --save --no-init-file < $rundir/jobs/merge.R\n";
    print MT "module unload R/3.2.2\n";
    close MT;
    my $hold_line = join ',',@hold_ids;
    if ( $hold_line eq '' ){
	print QSUB "\n##Merge HTSeq-count tables\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $mergeTables_job_id $rundir/jobs/$mergeTables_job_id.sh\n\n";
    } else {
	print QSUB "\n##Merge HTSeq-count tables\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $mergeTables_job_id -hold_jid $hold_line $rundir/jobs/$mergeTables_job_id.sh\n\n";
    }
    push (@hold_ids, $mergeTables_job_id);
}

if ( $opt{normalize} eq "yes"){
    #normalize the merged count table
    my ($normTables_job_id) = ("normTables_".get_job_id());
    open NM,">$rundir/jobs/$normTables_job_id.sh";
    print NM "\#!/bin/bash\n\n";
    print NM "uname -n > $rundir/logs/$normTables_job_id.host\n";
    #print NM "export MODULEPATH=/hpc/local/CentOS7/cog_bioinf/modules:\${MODULEPATH}\n";
    print NM "module load R/3.2.2\n";
    normalizeRscript($rundir);
    print NM "time R --save --no-init-file < $rundir/jobs/normalize.R\n";
    print NM "module unload R/3.2.2\n";
    close NM;
    my $hold_line = join ',',@hold_ids;
    if ( $hold_line eq '' ){
	print QSUB "\n##Normalize counts table\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $normTables_job_id $rundir/jobs/$normTables_job_id.sh\n\n";
    } else {
	print QSUB "\n##Normalize counts table\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $normTables_job_id -hold_jid $hold_line $rundir/jobs/$normTables_job_id.sh\n\n";
    }
}

if( $opt{rpkm} eq "yes"){
    #calculate rpkm's from merged count table
    my ($rpkm_job_id) = ("rpkm_".get_job_id());
    open RP, ">$rundir/jobs/$rpkm_job_id.sh";
    print RP "\#!/bin/bash\n\n";
    print RP "uname -n > $rundir/logs/$rpkm_job_id.host\n";
    #print RP "export MODULEPATH=/hpc/local/CentOS7/cog_bioinf/modules:\${MODULEPATH}\n";
    print RP "module load R/3.2.2\n";
    rpkmRscript($rundir);
    print RP "time R --save --no-init-file --args $rundir < $rundir/jobs/rpkm.R\n";
    print RP "module unload R/3.2.2\n";
    close RP;
    my $hold_line = join ',',@hold_ids;
    if ( $hold_line eq '' ){
	print QSUB "\n##Calculate RPKMs\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $rpkm_job_id $rundir/jobs/$rpkm_job_id.sh\n\n";    
    } else {
	print QSUB "\n##Calculate RPKMs\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $rpkm_job_id -hold_jid $hold_line $rundir/jobs/$rpkm_job_id.sh\n\n";    
    }
}

if ( $opt{de} eq "yes"){
    #perform differential expression analysis
    my ($de_job_id) = ("de_".get_job_id());
    open DE, ">$rundir/jobs/$de_job_id.sh";
    print DE "\#!/bin/bash\n\n";
    print DE "uname -n > $rundir/logs/$de_job_id.host\n";
    #print DE "export MODULEPATH=/hpc/local/CentOS7/cog_bioinf/modules:\${MODULEPATH}\n";
    print DE "module load R/3.2.2\n";
    deRscript($rundir);
    print DE "time R --save --no-init-file < $rundir/jobs/de.R\n";
    print DE "module unload R/3.2.2\n";
    close DE;
    my $hold_line = join ',',@hold_ids;
    if ( $hold_line eq '' ){
	print QSUB "\n##Perform DE analysis\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $de_job_id $rundir/jobs/$de_job_id.sh\n\n";    
    } else {
	print QSUB "\n##Perform DE analysis\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $de_job_id -hold_jid $hold_line $rundir/jobs/$de_job_id.sh\n\n";    
    }
}

if ( $opt{bamqc} eq "yes"){
    #invoke bamMetrics
    my ($bamQC_job_id) = ("bamQC_".get_job_id());
    push (@hold_ids, $bamQC_job_id);
    open BQ, ">$rundir/jobs/$bamQC_job_id.sh";
    print BQ "\#!/bin/bash\n\n";
    print BQ "uname -n > $rundir/logs/$bamQC_job_id.host\n";

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
    
    if ( $paired==0 ){
	print BQ "perl $opt{bamstats_path} -bam ".join(" -bam ",@final_bams)." -queue_threads 2 -debug -rna -ref_flat $opt{refflat_file} -ribosomal_intervals $opt{intervallist} -strand $picard_strand -single_end -genome $opt{fasta} -run_name $runname -output_dir $rundir/bamMetrics\n\n";
    }
    elsif ( $paired==1 ){
	print BQ "perl $opt{bamstats_path} -bam ".join(" -bam ",@final_bams)." -queue_threads 2 -debug -rna -ref_flat $opt{refflat_file} -ribosomal_intervals $opt{intervallist} -strand $picard_strand -genome $opt{fasta} -run_name $runname -output_dir $rundir/bamMetrics\n\n";
    }
    close BQ;
    if ( $opt{mapping} eq "no" ){
	print QSUB "\n##bamQC\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $bamQC_job_id $rundir/jobs/$bamQC_job_id.sh\n\n";
    }else{
	my $hold_line = join ',',@hold_mapping_ids;
	print QSUB "\n##bamQC\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $bamQC_job_id -hold_jid $hold_line $rundir/jobs/$bamQC_job_id.sh\n\n";
    }
    
    #cleaning
    print CL "if [ -f $rundir/bamMetrics/$runname.bamMetrics.pdf ]; then rm -r $rundir/bamMetrics/tmp; fi\n";
    close CL;
}

close SETTINGS;

my $hold_line = join ',',@hold_ids;
print QSUB "\n##cleaning\nqsub -l h_rt=02:00:00 -o $rundir/logs -e $rundir/logs -R yes -N $clean_job_id -hold_jid $hold_line $rundir/jobs/$clean_job_id.sh\n\n";
close QSUB;

system "sh $mainJobID" unless $opt{debug};



## SUBROUTINES ##
sub checkDesign{
    my $designfile = shift;
    my $status = 'good';
    my @array;
    open DESIGN, "<", $designfile or die "cannot open designfile\n";
    while ( my $line = <DESIGN> ){
	chomp($line);
	my @elements = split("\t",$line);
	$status = 'wrong' if $#elements+1 ne 2;
	push @array, $elements[1];
    }
    my @uniques = uniq @array;
    if ( $#uniques+1 ne 2 ){
	$status = 'wrong';
    }
    
    return($status);
}
sub mergeRscript{
    my $rundir = shift;
    my $runname = basename($rundir);
    open mergeRscript, ">$rundir/jobs/merge.R" or die "cannot open merge.R\n";
    print mergeRscript <<EOS;
dirs <- list.dirs(path="$rundir",recursive=F,full.names=F)
nr_cols=0
sampledir = ''
for ( dir in dirs ){
    sample=basename(dir)
    if ( (sample != 'read_counts') && (sample != 'bamMetrics') && (sample != 'logs') && (sample != 'jobs') && (sample != 'DEanalysis') && (sample != 'variantCalling') ){
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
    if ( (sample != 'read_counts') && (sample != 'bamMetrics') && (sample != 'jobs') && (sample != 'logs') && (sample != 'DEanalysis') && (sample != 'variantCalling') ) {
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
EOS
    close(mergeRscript);
}

sub normalizeRscript{
    my $rundir = shift;
    my $runname = basename($rundir);
    open normRscript, ">$rundir/jobs/normalize.R" or die "cannot open Rscript\n";
    print normRscript <<EOS;
library("DESeq")

outfile1="$rundir/read_counts/$runname\_readCounts_raw.txt"



#normalize raw read counts

counts=read.delim(outfile1, header=T, row.names=1)

conds <- factor(c(colnames(counts)))
cds <- newCountDataSet( counts, conds )
cds <- estimateSizeFactors( cds )

normalized_counts <- counts(cds, normalized=TRUE)

outfile2="$rundir/read_counts/$runname\_readCounts_normalized.txt"
write.table(normalized_counts, file=outfile2, row.names=T,col.names=NA, quote=F,sep="\\t")

EOS
    close(normRscript);
}

sub rpkmRscript{
    my $rundir = shift;
    my $runname = basename($rundir);
    my $gene_sizes = $opt{genesizes_file};
    open rpkmRscript, ">$rundir/jobs/rpkm.R" or die "cannot open Rscript\n";
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
write.table(df, file=outfile, row.names=T,col.names=NA, quote=F,sep="\t")

EOS
    close(rpkmRscript);
}

sub deRscript{
    my $rundir = shift;
    my $runname = basename($rundir);
    my $annotationdb = $annotationDB{$SPECIES};
    my $shortname = join('.', (split('\.',$annotationdb))[0,1]);
    my $design = abs_path($opt{design});
    open deRscript, ">$rundir/jobs/de.R" or die "cannot open Rscript\n";
    print deRscript <<EOS;
library($annotationdb)
library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)

outdir="$rundir/DEanalysis/"

# # ----- GENE ANNOTATIONS --------
# convert SYMBOL map to table
s = toTable($annotationdb.egSYMBOL)
# convert ENSEMBL map to table
t = toTable($annotationdb.egENSEMBL)
# merge symbol map table to ensembl map table
gene_annotations = merge(t, s, by.x = "gene_id", by.y = "gene_id")


# ----- IMPORT DATA --------
raw_read_counts <- read.table("$rundir/read_counts/$runname\_readCounts_raw.txt", sep="\\t",header=T,row.names=1,check.names=F)
orig.cols <- colnames(raw_read_counts)
colnames(raw_read_counts) <- make.names(colnames(raw_read_counts))
design = read.table("$design",row.names=1,header=F)
colnames(design) <- c('condition')
#put samples in same order as in raw_read_counts
conditions <- design[orig.cols,,drop=FALSE]


#---------- Standard design ----------

#construct a DESeqDataSet:
dds <- DESeqDataSetFromMatrix(countData = raw_read_counts, colData = conditions, design = ~ condition)

#Differential expression analysis
dds <- DESeq(dds)
dds <- dds[rowSums(counts(dds)) > 1,] #pre-filter low count genes
res <- results(dds) #access the results
resOrdered <- res[order(res\$padj),] #order by adjusted pvalue

#create MA plot: shows the log2 fold changes attributable to a given variable over the mean of normalized counts
#points are colored red if the adjusted p values is less than 0.1
#points which fall out of the window are plotted as open triangles pointing either up or down
png(paste(outdir,"$runname\_MAplot.png",sep=""))
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()


#---------- transformations ----------

# In order to test for differential expression, we operate on raw counts and use discrete distributions.
# However for other downstream analyses ( e.g.  for visualization or clustering ) it
# might be useful to work with transformed versions of the count data.

# The function rlog, stands for regularized log, transforming the original count data to the log2 scale
# by fitting a model with a term for each sample and a prior distribution on the coefficients which is
# estimated from the data.
rld <- rlog(dds)
rldMat <- assay(rld)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

#heatmap of the sample-to-sample distances
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hc <- hclust(distsRL)
png(paste(outdir,"$runname\_sampletosample_distances.png",sep=""))
heatmap.2(mat,Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col=rev(hmcol), margin=c(13,13))
dev.off()

#principal component plot of the samples
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
png(paste(outdir,"$runname\_PCAplot.png",sep=""))
ggplot(data, aes(PC1,PC2,color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()


# add gene symbols to results table
resOrdered\$id = rownames(resOrdered)
res1 = as.data.frame(resOrdered)
res2 = merge(res1, gene_annotations, by.x = "id", by.y = "ensembl_id", all.x=TRUE)
res4 = res2[order(res2\$pval),]
res5 = res4[which(!duplicated(res4\$id)),]
tablename <- paste(outdir,"$runname\_DEanalysis_all.txt",sep="")
write.table(res5, file=tablename, sep = "\\t", col.names = T, row.names = F, quote = F)


EOS
    close(deRscript);
}

sub get_job_id {
   my $id = tmpnam(); 
      $id=~s/\/tmp\/file//;
   return $id;
}

sub timestamp{
    my $t = localtime;
    return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",$t->year + 1900, $t->mon + 1, $t->mday, $t->hour, $t->min, $t->sec );
}
