library(DESeq2)

args<-commandArgs(TRUE)

mydata <- read.table(args[1], sep='\t', header = T, row.names = 1, stringsAsFactors = F)

print(mydata)

mock_data <- makeExampleDESeqDataSet(m=20)
colnames(mock_data) = rownames(mydata)

mtx = counts(mock_data)
mtx = rbind(mtx, 'mygene'=mydata$mygene)

annotations = data.frame(species=factor(mydata$species), diet = factor(mydata$diet))
rownames(annotations) = rownames(mydata)

dds = DESeqDataSetFromMatrix(countData=mtx, colData=annotations, design=~species+diet + species:diet)
dds = DESeq(dds)

print(colData(dds))
print(resultsNames(dds))

print('Model coefficients, log2 scale:')
print(coef(dds)['mygene',])

print('Species effect:')
r1 = results(dds, contrast=c('species', 'sp2', 'sp1'))
print(r1['mygene',])

print('Diet effect:')
r2 = results(dds, contrast=c('diet', 'A', 'B'))
print(r2['mygene',])

print('Interaction:')
r3 = results(dds, contrast=c(0,0,0,1))
print(r3['mygene',])

write.table(counts(dds), 'mocked_deseq2_data_with_simulated_gene.tsv', sep='\t', quote=F)
