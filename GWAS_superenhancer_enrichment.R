#last modified 04/26/2017
setwd("/Users/shuangliu/Box\ Sync/BrainspanNature")
options(stringsAsFactors=F)

#*********calculate PGC snp enrichment*************
snpscz=read.table('highpro_causalsnp.txt',header=F,sep='\t')#load schizophrenia snps for enrichment calculation
colnames(snpscz)=c('chr','start','end')

#*********define bedintersect function*************
bedTools.2in<-function(functionstring="bedIntersect",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
}

#************load  1kgp common snps as background*****************
load('snp1kgp.Rda')
kgpsubset=subset(snp1kgp,select=c('Chr','start','end'))
colnames(kgpsubset)=c('chr','start','end')
all_snp=rbind(snpscz,kgpsubset)
snp_list=all_snp[!duplicated(all_snp),]#remove GWAS snps from background

#************calculate enrichment p using fisher's exact test*******************
library("GenomicRanges") #http://davetang.org/muse/2013/01/02/iranges-and-genomicranges/
beddata=snpscz #test snps
colnames(beddata)=c('chr','start','end')
snpbed=with(beddata, GRanges(chr, IRanges(start, end)))
randset=snp_list#background snps
colnames(randset)=c('chr','start','end')
rand_snpbed=with(randset, GRanges(chr, IRanges(start, end)))
get_ens<-function(file_name){
  refdata=read.table(file_name,header=F)[,1:3];colnames(refdata)=c('chr','start','end')
  refdata=refdata[nchar(refdata$chr)<6,]
  refbed=with(refdata, GRanges(chr, IRanges(start, end)))
  snpsczF_overlap=sum(countOverlaps(snpbed,refbed))
  snpsczF_nooverlap=nrow(beddata)-snpsczF_overlap
  rand_overlap=sum(countOverlaps(rand_snpbed,refbed))
  rand_nooverlap=nrow(randset)-rand_overlap
  fishertable=cbind(rbind(snpsczF_overlap,snpsczF_nooverlap),rbind(rand_overlap,rand_nooverlap))
  colnames(fishertable)=c('causal','rand')
  rownames(fishertable)=c('mapped','not_mapped')
  ens_stats=fisher.test(fishertable,alternative ="greater")#chisq.test(fishertable)
  return(ens_stats$p.value)}

#******load H3K27ac peaks/enhancers/super enhancers or methylation sites*******
peaklist=c('PSNT.Adult','PSNT.infant','PSNT.Fetal','PSNT.Embryo','PSNT.DFC','PSNT.CBC')
peakname=paste(peaklist,'present-peaks.hg19.bedoverlapped.bed',sep='-')

enrichmentp=numeric()
for (i in 1:length(peakname))
{file_path=paste('/Users/shuangliu/Box\ Sync/BrainspanNature/overlap_enhancer/',peakname[i],sep='')
 enrichmentp[i]=get_ens(file_path)
}

bs_ens=data.frame(tissue=c(peaklist1,peaklist2),pval=enrichmentp)
save(bs_ens,file='overlap_enhancer/scz_enhancer_enrichment.Rda')
