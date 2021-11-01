#!/usr/bin/env Rscript

'path2genes prints out a list of genes starting from a pathway id (KEGG or Reactome).

Usage:
  path2genes (gene | reactome | kegg) <id> <out> [--affy=<a>]
  path2genes (-h | --help)
  path2genes (-v | --version)

Options:
  -h --help         Show this screen.
  -v --version      Show version.
  -a --affy=<a>     Use a table for AffyMetrix id conversion.
' -> doc

# Load the libraries
suppressPackageStartupMessages({
    library(docopt)
    library(data.table)
    #library(stringr)
    library(biomaRt)
})

args <- docopt(doc, version = 'degs2venn 0.1.0\n')
print(args)
source("algos.R")


# Actual workflow
write("Contacting the BioMart server...", stdout())
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")

if (args$gene == TRUE) {
    write("Running in Gene mode!", stdout())
    doGene(args, ensembl, out)
} else if (args$reactome == TRUE) {
    write("Running in Path mode with Reactome db!", stdout())
    dt <- doReactome(args, ensembl, out)
} else if (args$kegg == TRUE) {
    write("Running in Path mode with KEGG db!", stdout())
    dt <- doKEGG(args, ensembl, out)
}

if (length(args$affy) == 0) {
    write("Done!", stdout())
} else {
    write("Extracting the AffyMetix IDs!", stdout())
    doAffy(args, dt, out)
}

# Input/Output checks
#info <- test_input(args$csv)
#
#test_output(out)
#print(info)
