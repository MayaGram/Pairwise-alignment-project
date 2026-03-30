

def global_alignment(seq1:str, seq2: str)->int:

    match_score = 1
    mismatch_penalty = -1
    gap_penalty = -1

    prev = [j*gap_penalty for j in range(len(seq2)+1)]
    for i in range(1,len(seq1)+1):
        curr = [0]*(len(seq2)+1)

        curr[0] = i*gap_penalty

        for j in range(1,len(seq2)+1):
            match = prev[j-1] + (match_score if seq1[i-1]==seq2[j-1] else mismatch_penalty)
            delete = prev[j] +gap_penalty
            insert = curr[j-1]+gap_penalty

            curr[j] = max(match,delete,insert)

        prev = curr

    return prev[-1]




