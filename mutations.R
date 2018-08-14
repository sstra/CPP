

# NOTE: Will need to make functions specific to starting letter 
# p - protein, c - DNA, r - RNA, m - mitochondria, g - genome)

# It appears c, g, and m are all treated the same, in that they are all at the DNA level in some fashion


library(shiny)
library(DT) # requires development version for single row selection with datatables
library(DBI)
library(RMySQL)
library(ggplot2) # need development version for plotly (for horizontal bar)

#install_version("plotly", version = "4.6.0", repos = "http://cran.us.r-project.org")
library(plotly) # need development version 
library(stringr)



# NOTE: Need to add insertion-deletion function for all of the above

# Substitution: c|SUB|C|435|G -> c.C435G
# Deletion: p|DEL|508|F -> p.F508del
# Insertion: p|INS|1795|D -> p.1795insD
# Or
# Duplication: "c.1285-1301dup"	--> "c|DUP|1285_1301||	"c.1978(TATC)(1-2)"	--> "c|DUP|1978|TATC|1-2"  
# p.A3dup (single amino acid) or p.A3_S6dup (multiple amino acids)
# Frameshift: p|FS|G|46|| -> p.G46fsX
# Or: p|FS|A|31|V|41 -> p.A31VfsX41 (p.amino acid position amino acid fs terminates [x] at 41)

# nchar: To determine the length of a string
# "file:///C:/Users/SSTRA1275/Documents/R Tutorial Examples/CPP-master/mutation2pubtator.sample"
# f = read.delim("mutation2pubtator.sample", sep = "\t", header = TRUE)
f = read.delim("file:///C:/Users/SSTRA1275/Documents/R Tutorial Examples/CPP-master/mutation2pubtator.sample", sep = "\t", header = TRUE)
mutTerm <- as.character(f$Components)


sub <- "SUB"
del <- "DEL"
ins <- "INS"
dup <- "DUP"
fs <- "FS"
indel <- "INDEL"

#allMutations <- c("p|SUB|V|20|M", "p|SUB|R|433|W;RS#:56337013", "p|DEL|508|F", "|DEL||17")

#mutTerm <- allMutations
#orig <- mutTerm[3]

h <- grepl("^[[:alpha:]]\\|", mutTerm)
if (any(h)) {
  mutTerm <- mutTerm[h]
}

mutTerm <- gsub(";RS.+", "", mutTerm)

# Make RNA lowercase at end function
###########################################################################################################################################
# p|SUB|A|100|G p|SUB|A|100|

# Substitution format

formatSubstitutionAfterSplit <- function(x) {
  if (length(x) != 4) {
    stop("invalid format -- ", paste0(s, collapse = "|"))
  }
  return(paste0(x[1], x[3], x[2], ">", x[4]))
}

formatSubstitution <- function(sub, mutTerm, type) {
  if (type == "protein") {
    formatProteinSubstitution(sub, mutTerm) 
  } else if (type == "nucleotide") {
    formatDNASubsitution(sub, mutTerm)
  }
}


formatProteinSubstitution <- function(sub, mutTerm) {
  g <- grepl(sub, mutTerm)
  if (any(g)) {
    # look at first 'block' in p|SUB, etc
    # if sequence type is 'p', then :
    mutTerm[g] <- gsub("SUB", ".", mutTerm[g])
    mutTerm[g] <- gsub("\\|", "", mutTerm[g])
    # else if sequence type is 'c' then :
  }
  
  return(mutTerm)
}


formatNucleotideSubstitution <- function(sub, mutTerm) {
  g <- grepl(sub, mutTerm)
  if (any(g)) {
    mutTerm[g] <- gsub("\\|SUB", ".", mutTerm[g])
    s <- strsplit(mutTerm[g], "\\|")
    reformatted <- sapply(s, formatSubstitutionAfterSplit)
    mutTerm[g] <- reformatted 
  }
  
  return(mutTerm)
}


###########################################################################################################################################


# Deletion notationn

formatDeletionAfterSplit <- function(x) {
  if (length(x) != 4) {
    stop("invalid format -- ", paste0(s, collapse = "|"))
  }
  
  if (x[1] == "") {
    x[1] = "?"
  }
  return(paste0(x[1], ".", x[4], x[3], "del"))
} 


formatProteinDeletion <- function(del, mutTerm) {
  g <- grepl(del, mutTerm)
  if (any(g)) {
    s <- strsplit(mutTerm[g], "\\|")
    reformatted <- sapply(s, formatDeletionAfterSplit)
    mutTerm[g] <- reformatted 
  }
  
  return(mutTerm)
}


formatNucleotideDeletion <- function(del, mutTerm) {
  g <- grepl(del, mutTerm)
  if (any(g)) {
    mutTerm[g] <- gsub("\\|", "", mutTerm[g])
    mutTerm[g] <- gsub("DEL", ".", mutTerm[g])
  }
  
  return(mutTerm)
}


###########################################################################################################################################


# Insertion notation

formatProteinAndNucleoInsertion <- function(ins, mutTerm) {
  g <- grepl(ins, mutTerm)
  if (any(g)) {
    mutTerm[g] <- gsub("\\|INS\\|", ".", mutTerm[g])
    mutTerm[g] <- gsub("\\|", "ins", mutTerm[g])
  }
  
  return(mutTerm)
}


###########################################################################################################################################


# Deletion-Insertion notation

formatProteinAndNucleoDeletionInsertion <- function(indel, mutTerm) {
  g <- grepl(indel, mutTerm)
  if (any(g)) {
    mutTerm[g] <- gsub("\\|INDEL\\|", ".", mutTerm[g])
    mutTerm[g] <- gsub("\\|", "delins", mutTerm[g])
  }
  
  return(mutTerm)
}


###########################################################################################################################################


# Duplication notation

# Ambiguous amino acid: X; ambiguous base: N -> add addition specification to determine if protein, DNA, or RNA

formatProteinAndNucleoDuplication <- function(dup, mutTerm) {
  g <- grepl(dup, mutTerm)
  if (any(g)) {
  # mutTerm[g] <- gsub(";.+", "", mutTerm[g])
    mutTerm[g] <- gsub("\\|DUP\\|", ".", mutTerm[g])

    if(grepl("\\|\\|", mutTerm[g])) {
      m <- strsplit(mutTerm[g], "\\.")
      m <- unlist(m)
      if (m[1] == "p") {
        m[2] <- gsub("\\|\\|", "\\|X\\|", m[2])
        mutTerm[g] <- paste0(m[1], ".", m[2])
      }
      
      else {
        m[2] <- gsub("\\|\\|", "\\|N\\|", m[2])
        mutTerm[g] <- paste0(m[1], ".", m[2])
        
      }
    }
    
    if(grepl("_", mutTerm[g])) {
      mutTerm[g] <- gsub("_", "\\-", mutTerm[g])
      m <- strsplit(mutTerm[g], "\\|")
      m <- unlist(m)
      mutTerm[g] <- paste0(m[1], m[2], "dup")
    }
  
    if(grepl("[[:digit:]]$", mutTerm[g])) {
      m <- strsplit(mutTerm[g], "\\|")
      m < unlist(m)
      m[2] <- paste0("(", m[2], ")")
      m[3] <- paste0("(", m[3], ")")
      mutTerm[g] <- paste0(m[1], m[2], m[3])
    }
  }
  
  return(mutTerm)
}


###########################################################################################################################################


# Frameshift notation -> appears like this is only for protein when researching nomenclature

formatProteinFrameShift <- function(fs, mutTerm) {
  g <- grepl(fs, mutTerm)
  if(any(g)) {
  # mutTerm[g] <- gsub(";.+", "", mutTerm)
    mutTerm[g] <- gsub("\\|FS\\|", ".", mutTerm[g])
  
    if (grepl(mutTerm[g], "\\|\\|")) {
      mutTerm[g] <- gsub("\\|\\|", "fsX", mutTerm[g])
      mutTerm[g] <- gsub("\\|", "", mutTerm[g])
    }
  
    else {
      m <- strsplit(mutTerm[g], "\\|")
      m < unlist(m)
      m[length(m)] <- paste0("fsX", m[length(m)])
      mutTerm[g] <- paste0(m[1], m[2], m[3], m[4])
    }
  
  }
  
  return(mutTerm)
}


###########################################################################################################################################

#a <- cbind(orig, mutTerm)


# For application of the above functions:
# h finds every term with the specific label, then applies the corresponding functions to mutTerm[h]

h <- grepl("^p\\|", mutTerm)
if (any(h)) {
  # functions associated with protein
  
  mutTerm[h] <- formatProteinSubstitution(sub, mutTerm[h])
  mutTerm[h] <- formatProteinDeletion(del, mutTerm[h])
  mutTerm[h] <- formatProteinFrameShift(fs, mutTerm[h])
  
}


h <- grepl("^c\\||^g\\||^m\\||^r\\|", mutTerm)
if (any(h)) {
  # functions associated with DNA (genome and mitochondria seem to be synonomous)
  
  mutTerm[h] <- formatNucleotideSubstitution(sub, mutTerm[h])
  mutTerm[h] <- formatNucleotideDeletion(del, mutTerm[h])
  
}


h <- grepl("^[[:alpha:]]\\|", mutTerm)
if (any(h)) {
  mutTerm[h] <- formatProteinAndNucleoInsertion(ins, mutTerm[h])
  mutTerm[h] <- formatProteinAndNucleoDeletionInsertion(indel, mutTerm[h])
  mutTerm[h] <- formatProteinAndNucleoDuplication(dup, mutTerm[h])

}


h <- grepl("^r.", mutTerm)

mutTerm[h] <- tolower(mutTerm[h])

