#!/usr/bin/env Rscript
system("echo \n")
#### read in system arguements ######
options(echo=F) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# # testing purpose
#args<-c("test.fa",
#         "80",
#         "0.7",
#         "result",
#		"runtime",
#		"DIR")

# read in new SNP-free SARG ----

sargpath<-paste(args[6],"/ARGs_database_renamed.fnt.subset_subtype.name",sep="")

overlap.db<-read.delim(sargpath,stringsAsFactors = F,header=T)
colnames(overlap.db)<-c("subject","subtype",'type')

# read in taxonomy information of the marker genes
taxapath<-paste(args[6],"/taxa.info.RData",sep="")
load(taxapath)

if(length(args[4])==0) args[4]<-args[1] 



############
# ARG profile from NanoPore
# arg.coliform : LAST result before filtering
# arg.coliform2 : results filtered based on similarity > 70 and alignment > 50 length
# arg.colifom4: results after filtered out similar-position-ARG 
# arg.nanopore: results aggregated by ARG type
#############
library(plyr)
tmpname2<-paste(args[1],args[5],"sarg.tab",sep="_")
file<-paste("./tmp/",tmpname2,sep="")
arg.coliform<-read.delim(file,stringsAsFactors = F,header=F)

colnames(arg.coliform)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")
arg.coliform$query<-sapply(strsplit(arg.coliform$query,"-"),"[[",1)
# filtering hit based on similarity & alignment length of the ARG length
lookat<-which(arg.coliform$similarity>args[2] & arg.coliform$align.lenth/arg.coliform$s.len>args[3])
arg.coliform2<-arg.coliform[lookat,]

if(nrow(arg.coliform2)>0){
  #### remove the case where the same region on nanopore read hit to multiple ARGs ####
  tmp.lst<-split(arg.coliform2,arg.coliform2$query)
  
  # for one region exactly hit to multiple ARG, only the best hit (the one with highest bitscore also the first hit ) was kept
  tmp.lst<-lapply(tmp.lst, function(x) x[!duplicated(x$q.start),])
  tmp.lst<-lapply(tmp.lst, function(x) x[!duplicated(x$q.end),])
  # lapply(tmp.lst,nrow)
  
  # if one region hit to multiple ARG, then if the hited region is overlaped > 50% alignment length with the first hit (the hit with highest bitscore) then it will be removed, otherwise it will be kept.  ##
  
  for(g in 1:length(tmp.lst)) {
    x<-tmp.lst[[g]]
    # makesure q.start<q.end, flip q.end and q.start if not satisfy this standard
    lookat<-which(apply(x,1,function(y) as.numeric(y[7])> as.numeric(y[8])))
    tmp<-x[lookat,7]
    tmp2<-x[lookat,8]
    x[lookat,7]<-tmp2
    x[lookat,8]<-tmp
    # if nrow(x) > 2 then need to do clustering and then filter each cluster 
    if(nrow(x)>=2){
      x<-x[order(x$bitscore,decreasing=T),] # line with highest bitscore as first line
      tmp6<-list()
      tmp5<-vector() # store the line overlaped more than 80% with the first line, these lines should be deleted
      for(i in 1:(nrow(x)-1)){
        tmp5<-vector()
        for(j in (i+1):nrow(x)){
          tmp.start<-max(x[i,]$q.start,x[j,]$q.start)
          tmp.end<-min(x[i,]$q.end,x[j,]$q.end)
          overlap<-tmp.end-tmp.start
          # if no overlap with the first line， then overlap should be <=0, and this line should be kept for another loop
          # overlap with first line for more than 50% of alignment length
          if(overlap>x[j,]$align.lenth*0.5) {tmp5[j-1]<-j}
        }
        tmp6[[i]]<-tmp5
      }
      tmp6<-unique(unlist(tmp6))
      tmp6<-tmp6[!is.na(tmp6)]
      tmp.lst[[g]]<-x[-tmp6,]
    }
  }

  
  arg.colifom4<-rbind.fill(tmp.lst)
  
  
  # merge ARG annotation with ARDB type #
  tmp<-merge(arg.colifom4,overlap.db,by="subject")
  arg.colifom4<-tmp
  
}  else { 
  cat("Warning!: NO ARG identified\nbelow items will be empty；\narg.tab\narg.w.taxa.tab\n")
  arg.colifom4<-arg.coliform2
}


###############
# taxa profile
# taxa.info: taxaonomy information of marker genes
# taxa: taxa information of nanopore query after filter
###############
# read in the 2D.fa last marker gene result
tmpname3<-paste(args[1],args[5],"marker.tab",sep="_")
file2<-paste("./tmp/",tmpname3,sep="")

taxa<-read.delim(file2,stringsAsFactors = F,header=F)
colnames(taxa)<-c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len")
taxa$query<-sapply(strsplit(taxa$query,"-"),"[[",1)
# filtering hit based on similarity & alignment length of the marker gene length
lookat<-which(taxa$similarity> args[2] & taxa$align.lenth/taxa$s.len>args[3])
taxa<-taxa[lookat,]
taxa<-merge(taxa,taxa.info,by="subject")
taxa<-arrange(taxa,query,desc(bitscore))

# for each nanopore read only the taxa with the highest bitscore is kept。
taxa<-taxa[!duplicated(taxa$query),] # default deduplicated 就是取第一条

###############
# combine ARG profile with taxa profile
# arg.summary: the final results
################
# remove the unnecessary columns in taxa and arg.colifom4
taxa2<-taxa[,c("query","kindom","phylum","class","order","family","genus","species","sub.species")]

# get ARG profile of nanopore query with taxa classification
# results in arg.summary
if(nrow(arg.colifom4)>0){
	arg.w.taxa<-merge(arg.colifom4,taxa2,by="query")
	arg.w.taxa$species<-as.character(arg.w.taxa$species)
	arg.w.taxa<-arg.w.taxa[,c("query","subtype","type","kindom","phylum","class","order","family","genus","species","sub.species")]


	# --  write out ----
	# adjust colum sequence
	taxa<-taxa[,c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len","kindom","phylum","class","order","family","genus","species","sub.species")]
	arg.colifom4<-arg.colifom4[,c("query","subject","similarity","align.lenth","mismatch","gap","q.start","q.end","s.start","s.end","evalue","bitscore","s.len","q.len","subtype","type")]

	write.table(taxa,file=paste(args[4],"taxa.tab",sep="_"),
	quote=F,row.names = F,sep="\t")
	write.table(arg.colifom4,file=paste(args[4],"arg.tab",sep="_"),
	quote=F,row.names = F,sep="\t")
	write.table(arg.w.taxa,file=paste(args[4],"arg.w.taxa.tab",sep="_"),
	quote=F,row.names = F,sep="\t")

} else { 
  cat("Warning: NO ARG identified\nOnly taxa annotations were generated \n")
  write.table(taxa,file=paste(args[4],"taxa.tab",sep="_"),quote=F,row.names = F,sep="\t")
}

q(save="no")