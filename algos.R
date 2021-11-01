#!/usr/bin/env Rscript

doGene <- function(args, ensembl, out) {
    anno <- as.data.table(getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","reactome","kegg_enzyme","external_gene_name","description"), filters="ensembl_gene_id", values=args$id, mart=ensembl))
    fileName = file.path(args$out, paste(args$id, "gene2paths.tsv", sep="_"))
    write.table(anno, file=fileName, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
}

doReactome <- function(args, ensembl, out) {
    anno <- as.data.table(getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","reactome","kegg_enzyme","external_gene_name","description"), filters="reactome", values=args$id, mart=ensembl))
    fileName = file.path(args$out, paste(args$id, "Reactome", "path2genes.tsv", sep="_"))
    write.table(anno, file=fileName, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    return(anno)
}

doKEGG <- function(args, ensembl, out) {
    anno <- as.data.table(getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","reactome","kegg_enzyme","external_gene_name","description"), filters="kegg_enzyme", values=args$id, mart=ensembl))
    fileName = file.path(args$out, paste(args$id, "KEGG", "path2genes.tsv", sep="_"))
    write.table(anno, file=fileName, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    return(anno)
}

doAffy  <- function(args, dt, out) {
    clarioms <- fread(args$affy, skip=21)
    clarioms <- clarioms[,c(2,8)]
    #clarioms$gene_assignment <- str_split_fixed(clarioms$gene_assignment, "///", 5)[,3]
    #clarioms$affymetrix_symbol <- str_split_fixed(clarioms$gene_assignment, "//", 5)[,2]
    #clarioms$gene_assignment <- str_split_fixed(clarioms$gene_assignment, "//", 5)[,1]
    clarioms$gene_assignment <- sub(".*\\b(ENST\\w*).*", "\\1", clarioms$gene_assignment)
    setnames(clarioms, "gene_assignment", "ensembl_transcript_id")
    f <- merge(dt, clarioms, all.x=T)
    fileName = file.path(args$out, paste(args$id, "db", "path2genes.tsv", sep="_"))
    write.table(f, file=fileName, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    #f2 <- unique(f[,c(1,5,7)])
    f2 <- unique(f[,c(5,7)])
    fileName2 = file.path(args$out, paste(args$id, "unique", "db", "path2genes.tsv", sep="_"))
    write.table(f2, file=fileName2, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
}

# AnnotationDbi
#library(AnnotationDbi)
#library("clariomshumantranscriptcluster.db")
#keytypes(clariomshumantranscriptcluster.db)
#k <- head(keys(clariomshumantranscriptcluster.db, keytype="SYMBOL"))
#select(clariomshumantranscriptcluster.db, keys=k, columns=c("PROBEID","ENSEMBLTRANS","GENENAME"), keytype="SYMBOL")