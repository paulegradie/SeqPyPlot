args = commandArgs(trailingOnly=TRUE)

#args[1] = full genelist (from plot data)
#args[2] = output file name

plotterdata=read.csv(args[1], sep='\t')
genes=as.vector(plotterdata[2:length(plotterdata[,1]),1])

#Mouse
ensembl=useMart("ensembl",dataset="mmusculus_gene_ensembl")
conversionTable=getBM(attributes=c('entrezgene','mgi_symbol'),
                      filters='mgi_symbol',
                      values=genes,
                      mart=ensembl)

#Human
ensembl=useMart("ensembl",dataset="hsapeins_gene_ensembl")
conversionTable=getBM(attributes=c('entrezgene','hgnc_symbol'),
                      filters='mgi_symbol',
                      values=genes,
                      mart=ensembl)


write.table(conversionTable, file=args[2], sep = '\t', col.names=FALSE)