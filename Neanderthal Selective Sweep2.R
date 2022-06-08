library("rtracklayer")
library("ggplot2")

#Set working directory
setwd("~/PEC/Terry/Students/Asad/")

#Generate session for rtracklayer and select genome build
mySession <- browserSession ()
genome(mySession) <- "hg19"

#List available tracks on selected genome build
track.names <- trackNames(ucscTableQuery(mySession))
track.names

#Select appropriate data table from Neanderthal Selective Sweep Scan (S) track
track<-track.names[["Sel Swp Scan (S)"]]
# length(tableNames(ucscTableQuery(mySession, track=track))) #check tables in track
sweep<-getTable(ucscTableQuery(mySession,track=track,table=track))

# #Select appropriate data table from Ensembl genes track - NOT RUN
# track<-track.names[["Ensembl Genes"]]
# # tableNames(ucscTableQuery(mySession, track=track)) #check tables present on track
# genes<-getTable(ucscTableQuery(mySession,track=track,table="ensGene"))

#Select appropriate data table from Gencode track
track<-track.names[["GENCODE V39lift37"]]
tableNames(ucscTableQuery(mySession, track=track)) #check tables present on track
gencode_basic<-getTable(ucscTableQuery(mySession,track=track,table="wgEncodeGencodeBasicV39lift37"))
# gencode_attrs<-getTable(ucscTableQuery(mySession,track=track,table="wgEncodeGencodeAttrsV39lift37"))
# gencode_comp<-getTable(ucscTableQuery(mySession,track=track,table="wgEncodeGencodeCompV39lift37"))

#Generate table with appropriate columns from Gencode track data
gene_locations<-gencode_basic[,c("name","name2","cdsStart","cdsEnd")]

#Add 50k bases either side of each protein coding region
gene_locations$cdsStart_minus50k<-gene_locations$cdsStart-50000
gene_locations$cdsEnd_plus50k<-gene_locations$cdsEnd+50000

#Add columns of desired values, to be filled with data
gene_locations$length<-gene_locations$cdsEnd_plus50k-gene_locations$cdsStart_minus50k
gene_locations$sweep_count<-NA #Number of Neanderthal sweep hits within region
gene_locations$sweep_min<-NA #Minimum sweep value within region
gene_locations$sweep_med<-NA #Median sweep value
gene_locations$sweep_mean<-NA #Mean sweep value
gene_locations$sweep_max<-NA #Max sweep value
gene_locations$sweep_sd<-NA #SD of sweep values

#Iterate through Gencode track data and generate results for each region
for(i in 1:nrow(gene_locations)){
  v<-which(findInterval(sweep$start, c(gene_locations$cdsStart_minus50k[i], gene_locations$cdsEnd_plus50k[i]))==1)
  gene_locations$sweep_count[i]<-length(v)
  gene_locations$sweep_min[i]<-min(sweep$value[v]/gene_locations$length[i])
  gene_locations$sweep_med[i]<-median(sweep$value[v]/gene_locations$length[i])
  gene_locations$sweep_mean[i]<-mean(sweep$value[v]/gene_locations$length[i])
  gene_locations$sweep_max[i]<-max(sweep$value[v]/gene_locations$length[i])
  gene_locations$sweep_sd[i]<-sd(sweep$value[v]/gene_locations$length[i])
}

#Generate sample of gene regions
v<-sample(c(1:nrow(gene_locations)),replace = T,size = 10000)
res<-data.frame(gene=gene_locations$name[v],min_S=gene_locations$sweep_min[v],Dataset=rep("Genome Sample",10000))

#Plot distribution of minimum S value
p1<-ggplot(res,aes(y=min_S,x=Dataset))+
  geom_violin(size=1)+
  geom_boxplot(size=1,width=.5)+
  ylab(label = "Minimum Length-Normalized S Value")+
  xlab(label = "Dataset")+
  theme_bw(base_size=22)

#Save plot
png("Distribution_of_Genome_Sweep.png",height=800,width=800)
p1
dev.off()

#Save data
write.csv(res,"Genome_sampled_s_values.csv",row.names = T)
write.csv(gene_locations,"Full_genome_w_Sweep.csv",row.names = T)

#Growth Genes
growth_genes<-read.csv("Growth_genes_by_year.csv",header=T)

growth_gene_locations<-list()
for(i in 1:length(unique(growth_genes$Year))){
growth_gene_locations[[i]]<-gene_locations[na.omit(which(!is.na(match(gene_locations$name2,growth_genes$Gene[growth_genes$Year==unique(growth_genes$Year)[i]])))),]
}

for(j in 1:length(unique(growth_genes$Year))){
  locs<-growth_gene_locations[[j]]
  #Iterate through Gencode track data and generate results for each region
  for(i in 1:nrow(locs)){
    v<-which(findInterval(sweep$start, c(locs$cdsStart_minus50k[i], locs$cdsEnd_plus50k[i]))==1)
    locs$sweep_count[i]<-length(v)
    locs$sweep_min[i]<-min(sweep$value[v]/locs$length[i])
    locs$sweep_med[i]<-median(sweep$value[v]/locs$length[i])
    locs$sweep_mean[i]<-mean(sweep$value[v]/locs$length[i])
    locs$sweep_max[i]<-max(sweep$value[v]/locs$length[i])
    locs$sweep_sd[i]<-sd(sweep$value[v]/locs$length[i])
  }
  growth_gene_locations[[j]]<-locs
}

growth_gene_res<-do.call(rbind,growth_gene_locations)

growth_gene_res$Dataset<-rep(unique(growth_genes$Year),lapply(growth_gene_locations,nrow))
growth_gene_res<-data.frame(gene=growth_gene_res$name,min_S=growth_gene_res$sweep_min,Dataset=growth_gene_res$Dataset)

#Plot distribution of minimum S value
p1<-ggplot(growth_gene_res,aes(y=min_S,x=Dataset))+
  geom_violin(size=1)+
  geom_boxplot(size=1,width=.3)+
  ylab(label = "Minimum Length-Normalized S Value")+
  xlab(label = "Dataset")+
  theme_bw(base_size=22)

#Save plot
png("Growth_genes_by_year_sweep.png",height=800,width=800)
p1
dev.off()

#Save data
write.csv(growth_gene_res,"Growth_genes_by_year_s_values.csv",row.names = T)
