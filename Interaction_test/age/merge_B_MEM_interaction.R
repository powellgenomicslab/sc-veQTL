####
library(data.table)

m=c()
v=c()
d=c()

for(i in 1:22){
	tmpm=fread(paste0("./B_MEM_cell_5_mean_mx/OneK1K_B_MEM_chr",i,".cis_qtl_top_assoc.txt.gz"),header=T)
	tmpv=fread(paste0("./B_MEM_cell_5_variance_mx/OneK1K_B_MEM_chr",i,".cis_qtl_top_assoc.txt.gz"),header=T)
	tmpd=fread(paste0("./B_MEM_cell_5_dispersion_mx/OneK1K_B_MEM_chr",i,".cis_qtl_top_assoc.txt.gz"),header=T)

	m=rbind(m,tmpm)
	v=rbind(v,tmpv)
	d=rbind(d,tmpd)
}

m=as.data.frame(m)
v=as.data.frame(v)
d=as.data.frame(d)

write.table(m,"B_MEM_mean_interaction.txt",row.names=F,col.names=T,quote=F)
write.table(v,"B_MEM_variance_interaction.txt",row.names=F,col.names=T,quote=F)
write.table(d,"B_MEM_dispersion_interaction.txt",row.names=F,col.names=T,quote=F)

gene=c("HLA-DQA1","RPS26","SMDT1")
m2=m[m$phenotype_id %in% gene,]
v2=v[v$phenotype_id %in% gene,]
d2=d[v$phenotype_id %in% gene,]
res=rbind(m2,v2,d2)

write.table(res,"B_MEM_candidate_interaction.txt",row.names=F,col.names=T,quote=F)

m=read.table("B_IN_mean_interaction.txt",header=T)
v=read.table("B_IN_variance_interaction.txt",header=T)
d=read.table("B_IN_dispersion_interaction.txt",header=T)
gene=c("HLA-DQA1","RPS26","SMDT1")
m2=m[m$phenotype_id %in% gene,]
v2=v[v$phenotype_id %in% gene,]
d2=d[v$phenotype_id %in% gene,]

res=rbind(m2,v2,d2)
write.table(res,"B_IN_candidate_interaction.txt",row.names=F,col.names=T,quote=F)

####
