library('ape')
library('RColorBrewer')
library('rgl')
source('lib.r')

# returns index of root of tree
"root.of" <- function(tree) return(length(tree$tip.label) + 1)

cols <- c(brewer.pal(9,'Set1')[-6],brewer.pal(8,'Set2'),brewer.pal(9,'Set3')); cols <- sprintf('%saa',cols)

cat('Loading fasta\n')
gg <- read.FASTA('97_otus_p.00125.fasta')
# gg <- read.FASTA('97_otus_p.000125.fasta')

cat('Building distances\n')
dd <- as.matrix(dist.dna(gg))

cat('   making tree\n')
ggt <- bionj(dd)
ntips <- length(ggt$tip.label)

cat('   getting tree distances\n')
ddt <- dist.nodes(ggt)

# read in the taxonomy names and extract just the phylum
tax <- read.table('97_otu_taxonomy.txt',sep='\t')
tax <- tax[match(ggt$tip.label,tax[,1]),]
phylum <- sapply(strsplit(as.character(tax[,2]),'; '),function(xx) substr(xx[2],4,nchar(xx[2])));
phylum[phylum %in% names(table(phylum))[table(phylum) <= 5]] <- 'Other'


res <- embed.phylo(ggt,topo.weight=10,verbose=TRUE,optim.method='BFGS')
pc.sphere <- spherify.coords(ggt,res$coords,center.type='root',eps=10)
plot.phylo.3D(ggt, pc.sphere ,phylum)


# plot a fan version of the tree
pdf('tree.pdf',width=5, height=5)
plot(ggt,'fan',tip.color=cols[1:length(unique(phylum[1:ntips]))][as.numeric(as.factor(phylum[1:ntips]))], cex=.3)
dev.off()

# save STL
writeSTL('tree.stl',lineRadius=.003,pointRadius=.007)

# save webGL version
writeWebGL('webgl')

# save the tree
write.tree(ggt,file='16s_tree.tre')

