import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def count_matches(sub_seq1, sub_seq2):
    matches = 0

    for i in range(min(len(sub_seq1), len(sub_seq2))):
        if (sub_seq1[i] == sub_seq2[i]):
            matches += 1

    return matches

def slice_sequence(sequence, window_size, step):
    slices  = []
    start = 0
    while start + window_size <= len(sequence):
        slice = sequence[start:start + window_size]
        slices.append(slice)
        start += step
    return slices

def compare_sub_sequence(sub_seq1, sub_seq2, threshold):
    matches = count_matches(sub_seq1, sub_seq2)

    if (matches >= threshold):
        return 1
    else:
        return 0

def dot_matrix(seq1, seq2, window_size = 1, allowed_mismatches = 0, step = 1):

    seq1_sliced_sequences = slice_sequence(seq1, window_size, step)
    seq2_sliced_sequences = slice_sequence(seq2, window_size, step)

    rows = len(seq1)
    columns = len(seq2)

    matrix = np.zeros((rows, columns), dtype=int)

    for seq1_window, column_index in zip(seq1_sliced_sequences, range(len(seq1))):
        for seq2_window, row_index in zip(seq2_sliced_sequences, range(len(seq2))):
            row_index_for_dot = int(row_index * step + window_size / 2)
            column_index_for_dot = int(column_index * step + window_size / 2)

            matrix[column_index_for_dot, row_index_for_dot] = compare_sub_sequence(seq1_window, seq2_window, allowed_mismatches)
    return matrix

def plot_dot_matrix(matrix):
    cmap = LinearSegmentedColormap.from_list('', ['white', 'black'])

    graph, plt1 = plt.subplots(figsize = (100,100), dpi= 240)

    # plt1.hlines(y=range(len(seq1)),
    #       xmin=0, xmax=len(seq1),
    #       color='gray', alpha=0.7, linewidth=1,linestyles='dashdot')
    # plt1.vlines(x=range(len(seq2)),
    #       ymin=0, ymax=len(seq2),
    #       color='gray', alpha=0.7, linewidth=1,linestyles='dashdot')

    plt1.imshow(matrix, cmap=cmap)
    plt1.set_xticks(range(len(seq1)))
    plt1.set_yticks(range(len(seq2)))
    plt1.invert_yaxis()
    plt1.set_yticklabels(list(seq2))
    plt1.set_xticklabels(list(seq1))
    plt.show()


parser = argparse.ArgumentParser()

parser.add_argument(
    '-s',
    '--sequences',
    help = 'Enter the sequences to be compared',
    nargs = 2
)

parser.add_argument(
    '-w',
    '--window_size',
    help = 'Enter the Window size of dot-matrix plot',
    default = 1
)

parser.add_argument(
    '-step',
    '--step_size',
    help = 'Enter the Step size of dot-matrix plot',
    default = 1
)

parser.add_argument(
    '-t',
    '--threshold',
    help = 'Enter the allowed number of matches taken as threshold',
    default = 1
)

args = parser.parse_args()


if args.sequences:
    seq1 = args.sequences[0]
    seq2 = args.sequences[1]

    matrix = dot_matrix(seq1, seq2, int(args.window_size), int(args.threshold), int(args.step_size))

    plot_dot_matrix(matrix)

