import numpy as np

def count_mutations(seq, wt_seq):
    """
    Counts the number of mutations in an amino acid sequence compared to a wild-type sequence
    using the Needleman-Wunsch algorithm to globally align the sequences.
    
    Args:
    seq (str): The amino acid sequence to compare.
    wt_seq (str): The wild-type amino acid sequence to compare against.
    
    Returns:
    int: The number of mutations in the sequence compared to the wild-type sequence.
    """
    # Create a scoring matrix
    match_score = 1
    mismatch_score = -1
    gap_score = -1
    nrow, ncol = len(seq) + 1, len(wt_seq) + 1
    score_matrix = np.zeros((nrow, ncol), dtype=int)
    score_matrix[:, 0] = gap_score * np.arange(nrow)
    score_matrix[0, :] = gap_score * np.arange(ncol)
    for i in range(1, nrow):
        for j in range(1, ncol):
            match = score_matrix[i-1, j-1] + (match_score if seq[i-1] == wt_seq[j-1] else mismatch_score)
            delete = score_matrix[i-1, j] + gap_score
            insert = score_matrix[i, j-1] + gap_score
            score_matrix[i, j] = max(match, delete, insert)
    
    # Trace back the optimal alignment
    i, j = nrow-1, ncol-1
    mutations = 0
    while i > 0 or j > 0:
        if i > 0 and score_matrix[i, j] == score_matrix[i-1, j] + gap_score:
            i -= 1
            mutations += 1
        elif j > 0 and score_matrix[i, j] == score_matrix[i, j-1] + gap_score:
            j -= 1
            mutations += 1
        else:
            i -= 1
            j -= 1
            if seq[i] != wt_seq[j]:
                mutations += 1
    
    return mutations
