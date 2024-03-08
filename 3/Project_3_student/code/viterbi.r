#!/usr/bin/env Rscript

#' @param E n times m matrix with the emission log probabilities of the n latent variables and m observed variables
#' @param Tr n times n matrix with the transition log probabilities
#' @param I numeric vector of length n with the initial log probabilities of the n latent variables
#' @param p a data.frame with a Variable named AminoAcids which holds the amino acid sequence as a character string
viterbi <- function(E, Tr, I, p) {
    .as.array <- function(.) stringr::str_split(., "")[[1]]
    unique.ss <- c("B", "C", "E", "G", "H", "I", "S", "T")
    unique.aa <- c("A", "C", "D", "E", "F", "G", "H", "I",
                   "K", "L", "M", "N", "P", "Q", "R", "S",
                   "T", "U", "V", "W", "X", "Y")
    for (k in seq(nrow(p))) {
        sequence <- p$AminoAcids[k]
        aa.vec <- .as.array(sequence) %>% match(unique.aa)
        P      <- matrix(0, nrow(E), length(aa.vec))
        Ptr    <- matrix(0, nrow(E), length(aa.vec))
        
        ## sets the paths
        for (i in seq(length(aa.vec))) {
            if (i == 1) {
                P[, i] <- I + E[, aa.vec[i]]
            } else {
                for (j in seq(nrow(E))) {
                    p.loc    <- P[, i - 1] + Tr[, j] + E[j, aa.vec[i]]
                    P[j, i] <- max(p.loc)
                    Ptr[j, i] <- which.max(p.loc)
                }
            }
        }
        
        ## backtrace: computes the most likely path
        Phi <- vector(mode="integer",   length=length(aa.vec))
        Phi[length(Phi)] <- which.max(P[, ncol(P)])
        ## we start at the back, just as with Needleman-Wunsch or Smith-Waterman
        for (i in seq(from=length(aa.vec), to=2)) {
            Phi[i - 1] <- Ptr[Phi[i], i]
        }
        
        states <- unique.ss[Phi]
        p$PredictedStructure[k] <- paste(states, collapse="")
    }
    return(p)
}

assign_row_names <- function(M) {
    row.names(M) <- c("C", "S", "B", "E", "T", "I", "H", "G")
    return(M)
}

assign_col_names <- function(M, aa=TRUE) {
    if (aa) {
        colnames(M) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y")
    } else {
        colnames(M) <- c("C", "S", "B", "E", "T", "I", "H", "G")
    }
    
    return(M)
}

init_matrices <- function(data) {
    # Initialisation of I, T, and E matrices using maximum likelihood

    # I - 8 x 1 matrix - initial state probabilities
    I <- matrix(0, nrow=8, ncol=1)
    # T - 8 x 8 matrix - transition probabilities
    T <- matrix(0, nrow=8, ncol=8)
    # E - 8 x 20 matrix - emission probabilities
    E <- matrix(0, nrow=8, ncol=22)

    I <- assign_row_names(I)
    T <- assign_row_names(T)
    E <- assign_row_names(E)

    T <- assign_col_names(T, aa=FALSE)
    E <- assign_col_names(E)

    I <- initial_states(data, I)
    print(E)
}

initial_states <- function(data, M) {
    # Get the third column of the dataframe
    column <- data[, 3]
    
    # Initialize a frequency table
    frequency <- table(substr(column, 1, 1))
    
    # Write values from frequency table into I matrix
    M[row.names(M) %in% names(frequency), 1] <- as.numeric(frequency)

    # Return initial states matrix
    return(M)
}

#' Main function that accepts command-line arguments
#' @param args a character vector of command-line arguments
main <- function(args) {
  # Parse command-line arguments
  if (length(args) != 1) {
    stop("Invalid number of arguments. Usage: Rscript viterbi.r <data_folder>")
  }
  data_folder <- args[1]
  
  # Loading data
  prot_train_df <- read.csv(paste0(data_folder, "proteins_train.tsv"), header=FALSE, sep="\t")
  prot_test_df <- read.csv(paste0(data_folder, "proteins_test.tsv"), header=FALSE, sep="\t")
  prot_new_df <- read.csv(paste0(data_folder, "proteins_new.tsv"), header=FALSE, sep="\t")

  init_matrices(prot_train_df)

}

args = commandArgs(trailingOnly=TRUE)
main(args)

