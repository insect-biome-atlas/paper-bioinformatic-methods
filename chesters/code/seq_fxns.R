library(ape)
library(seqinr)

# Function for renaming sequences with tip labels according to
# metadata file. Matching is based on GenBank record ID
rename_seqs <- function(seqs, meta) {

    x <- names(seqs)
    w <- character()
    
    for (i in 1:length(x)) {
        w[i] <- strsplit(x[i]," ")[[1]][1]
    }

    w <- unlist(w)

    names(seqs) <- meta$TipLabel[match(w,meta$NCBI)]

    seqs
}

# Function for extracting coding sequences in DNAbin format (ape)
extract_coding <- function(seqs, meta) {

    x <- names(seqs)

    stop_codon1 <- c("18","88","88")   # TAA in DNAbin
    stop_codon2 <- c("18","88","28")   # TAG in DNAbin

    for (i in 1:length(x)) {
        
        meta_idx <- match(x[i],meta$TipLabel,nomatch=NA)
#        cat ("i=", i, "; meta_idx=", meta_idx, "\n")
        if (is.na(meta_idx))
            cat("No match for ",x[i],"\n")
        
        len <- meta$Stop[meta_idx] - meta$Start[meta_idx] + 1
        if (len %% 3 != 0)
            cat ("Length for sequence", i, "not divisible by 3\n")
        last_codon <- as.character(seqs[[i]][(meta$Stop[meta_idx]-2):meta$Stop[meta_idx]])
        if (identical(last_codon,stop_codon1) || identical(last_codon,stop_codon2))
            seqs[[i]] <- seqs[[i]][meta$Start[meta_idx]:(meta$Stop[meta_idx]-3)]
        else
            seqs[[i]] <- seqs[[i]][meta$Start[meta_idx]:meta$Stop[meta_idx]]

    }

    seqs
}

# Function for converting aligned fasta sequences to nexus data block
fasta2nexus <- function(infile, outfile) {

    seqs <- read.alignment(infile, format="fasta")

    ntax <- length(seqs$seq)
    nchar <- nchar(seqs$seq[1])

    cat ("#NEXUS\n\nbegin data;\n\tdimensions ntax =", ntax, " nchar =", nchar, ";\n", sep="", file=outfile)
    cat ("\tformat datatype=dna gap=- interleave=no;\n\tmatrix\n", sep="", file=outfile, append=TRUE)

    for (i in 1:length(seqs$seq)) {
        cat(seqs$nam[i], "\t", toupper(seqs$seq[i]), "\n", sep="", file=outfile, append=TRUE)
    }

    cat("\t;\nend;\n", file=outfile, append=TRUE)
}

# Function for finding stop codons in aligned DNA sequences
# in ape format, using trans fxn in ape
no_stop_codons <- function(seqs, codonstart) {

    aa_seqs <- trans(seqs,code=5, codonstart=codonstart)
    x <- logical()
    for (i in 1:length(aa_seqs))
        x[i] <- is.na(match("*",as.character(aa_seqs[[i]])))

    x
}

