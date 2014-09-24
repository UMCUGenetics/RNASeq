#Run by typing: R --vanilla --args rundir < merge_count_tables.r

#Author: Annelies Smouter
#Date: may 20, 2014

args <- commandArgs(TRUE)
rundir <- args[1]

library("DESeq")

dirs <- list.dirs(path=rundir,recursive=F,full.names=F)
nr_cols=0
sampledir = ''
for ( dir in dirs ){
    sample=basename(dir)
    if ( sample != 'read_counts' ){
	nr_cols = nr_cols+1
	sampledir = dir
    }
}

nr_rows <- nrow(read.table(paste(sampledir,'/read_counts/',basename(sampledir),'_htseq_counts.txt',sep=""),sep="\t",header=F))-5

output <- matrix(ncol=nr_cols+1, nrow=nr_rows)
col_count=1
samplenames=c('gene')

for ( dir in dirs ) {
    sample=basename(dir)
    if ( sample != 'read_counts' ) {
	samplenames <- append(samplenames, as.character(sample))
	countsfile=paste(dir,'/read_counts/',sample,'_htseq_counts.txt',sep="")
	counts=read.table(countsfile,sep="\t",header=F)
	colnames(counts) <- c('gene',sample)
	genes<-counts$gene

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

outfile1=paste(rundir,'/read_counts/raw_read_counts_table.txt', sep="/")
write.table(merged_table,file=outfile1,row.names=F,col.names=T,quote=F,sep="\t")



#normalize raw read counts

counts=read.delim(outfile1, header=T, row.names=1)

conds <- factor(c(colnames(counts)))
cds <- newCountDataSet( counts, conds )
cds <- estimateSizeFactors( cds )

normalized_counts <- round(counts(cds, normalized=TRUE))

outfile2= paste(rundir,"read_counts/normalized_read_counts_table.txt",sep="/")
write.table(normalized_counts, file=outfile2, row.names=T,col.names=NA, quot=F,sep="\t")
