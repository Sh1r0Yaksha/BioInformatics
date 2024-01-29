rm(list = ls())
# Install these libraries in case they don't exist
library("stringr")
library("stringi")

# Put Sequences here
seq1 <- "ATTTCTAGGGAAGTTTCAGATGGAGAAATGTGTGGACCATTCAGGCCCTAAAAGTTGGTATAAAGAT"
seq2 <- "AAGCACCGCGTAGGCCAGCTGGCCGGATCCCGCCGTCTGTCATGGCGGCCCCCATCCTGAAAGG"

# The respective match score, mismatch score and gap penalty
match_score <- 1
mismatch_score <- -1
gap_penalty <- -2

# Breaking the string into a vector of characters - like ['A','C','T','G']
s1_v <- strsplit(seq1, NULL)[[1]]
s2_v <- strsplit(seq2, NULL)[[1]]

# Function to check whether two base pairs match or don't, and returns the respective score
check_match <- function(a, b, match, mismatch) {
  if (str_trim(a) == str_trim(b)) {
    return(match)
  } else {
    return(mismatch)
  }
}

# Function to generate the similarity matrix based on gap penalty and match/mismatch score
generate_matrix <- function(seq1, seq2, match_score, mismatch_score, gap_score) {

  sim_matrix <- matrix(0, nrow = length(seq2) + 1, ncol = length(seq1) + 1)

  for (i in 1 : ncol(sim_matrix)) {
    sim_matrix[1, i] <- gap_score * (i - 1)
  }
  for (i in 1: nrow(sim_matrix)) {
    sim_matrix[i, 1] <- gap_score * (i - 1)
  }
  for (seq1_col in 1:length(seq1)) {
    for (seq2_row in 1:length(seq2)) {
      match <- sim_matrix[seq2_row, seq1_col] +
        check_match(seq2[seq2_row], seq1[seq1_col],
                    match_score, mismatch_score)
      delete <- sim_matrix[seq2_row, seq1_col + 1] + gap_score
      insert <- sim_matrix[seq2_row + 1, seq1_col] + gap_score
      sim_matrix[seq2_row + 1, seq1_col + 1] <- max(match, delete, insert)
      rm(match, insert, delete)
    }
  }
  return(list(sim_matrix = sim_matrix,
              score = sim_matrix[nrow(sim_matrix), ncol(sim_matrix)]))
}

# Function to align the two sequences tracing back the similarity matrix
pairwise_align <- function(seq1, seq2, match_score, mismatch_score, gap_penalty) {

  result <- generate_matrix(seq1, seq2, match_score, mismatch_score, gap_penalty)
  sim_matrix <- result$sim_matrix
  alignment_a <- ""
  alignment_b <- ""
  col_seq1 <- length(seq1) + 1
  row_seq2 <- length(seq2) + 1

  while (col_seq1 > 1 || row_seq2 > 1) {
    if (col_seq1 > 1 && row_seq2 > 1) {
      diag <- sim_matrix[row_seq2, col_seq1] == sim_matrix[row_seq2 - 1,col_seq1 - 1] +
              check_match(seq2[row_seq2 - 1], seq1[col_seq1 - 1], match_score, mismatch_score)
      delete <- sim_matrix[row_seq2, col_seq1] == sim_matrix[row_seq2 - 1, col_seq1] + gap_penalty
      insert = sim_matrix[row_seq2, col_seq1] == sim_matrix[row_seq2, col_seq1 - 1] + gap_penalty

      if (insert) {
        alignment_a <- paste0(alignment_a, seq1[col_seq1 - 1])
        alignment_b <- paste0(alignment_b, "-")
        col_seq1 <- col_seq1 - 1
      } else if (delete) {
        alignment_a <- paste0(alignment_a, "-")
        alignment_b <- paste0(alignment_b, seq2[row_seq2 - 1])
        row_seq2 <- row_seq2 - 1
      } else if (diag) {
        alignment_a <- paste0(alignment_a, seq1[col_seq1 - 1])
        alignment_b <- paste0(alignment_b, seq2[row_seq2 - 1])
        col_seq1 <- col_seq1 - 1
        row_seq2 <- row_seq2 - 1
      } else {
        return()
      }
    } else if (col_seq1 > 0) {
      alignment_a <- paste0(alignment_a, seq1[col_seq1 - 1])
      alignment_b <- paste0(alignment_b, "-")
      col_seq1 <- col_seq1 - 1
    } else if (row_seq2 > 0) {
      alignment_a <- paste0(alignment_a, "-")
      alignment_b <- paste0(alignment_b, seq2[row_seq2 - 1])
      row_seq2 <- row_seq2 - 1
    } else {
      return()
    }
  }

  alignment_a <- stri_reverse(alignment_a)
  alignment_b <- stri_reverse(alignment_b)


  return(list(alignment_a = alignment_a, alignment_b = alignment_b, score = result$score))
}

# Printing the results
align <- pairwise_align(s1_v, s2_v, match_score, mismatch_score, gap_penalty)

cat("Alignment 1: ", align$alignment_a, "\n")
cat("Alignment 2: ", align$alignment_b, "\n")
cat("Score:", align$score, "\n")