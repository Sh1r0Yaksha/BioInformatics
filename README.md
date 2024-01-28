# BioInformatics

Few BioInformatics Algorithm written in Python and R

The algorithms are - 

* Dot Matrix plot for sequence alignment
* Needleman-Wunsch algorithm for global alignment of two sequences

# Files

## /Python

### dot_matrix.ipynb

A jupyter notebook for plotting the dot_matrix alignment of two sequences with given window size, steps and threshold

Enter the sequences in the variables - `seq1` and `seq2` and the window size, steps and threshold in the variables `window_size`, `step` and `threshold` respectively and run all cells to get the plot

### dot_matrix.py

Python file which you can run and provide commands to get the dot matrix plot, same as the jupyter notebook with following arguments -

* **-s** or **--sequences**: The two sequences to be compared
* **-w** or **--window_size**: Window size of the alignment
* **-step** or **--step_size**: Step size of the plot
* **-t** or **--threshold**: Allowed number of matches in a window

Example: `python3 dot_matrix.py -s ACTGTGGACT GTAGTCAGT -w 3 -step 2 -t 2`

### needleman_wunsch.ipynb

A jupyter notebook for finding the global sequence alignment of two sequences with given match score, mismatch score and gap penalty

Enter the sequences in the variables - `seq1` and `seq2` and the match score, mismatch score and gap penalty in the variables `match_score`, `mismatch_score` and `gap_penalty` respectively and run all cells to get the alignment

### needleman_wunsch.py

Python file which you can run and provide commands to get the global sequence alignment, same as the jupyter notebook with following arguments -

* **-s** or **--sequences**: The two sequences to be compared
* **-match** or **--match_score**: Score given when two base pairs match
* **-mismatch** or **--mismatch_score**: Score given when two base pairs do not match
* **-g** or **--gap_penalty**: Score given when there is a gap introduced in the two alignments

Example: `python3 needleman_wunsch.py -s ACTGTGGACT GTAGTCAGT -match 1 -mismatch -1 -g -2`

## /R

### needleman_wunsch.r

R file for finding the global sequence alignment of two sequences with given match score, mismatch score and gap penalty

Enter the sequences in the variables - `seq1` and `seq2` and the match score, mismatch score and gap penalty in the variables `match_score`, `mismatch_score` and `gap_penalty` respectively and run the file to get the alignment

# Note

* Some of the code might be inefficient or even wrong.
* The `needleman_wunsch.r` file runs very slow and gives unaccurate results if the sequences provided are very large (more than 60 base pairs)
