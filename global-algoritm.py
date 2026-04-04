

def nw_score(seq1, seq2):
    gap_penalty = -1
    match_score = 1
    mismatch_penalty = -1

    prev = [j*gap_penalty for j in range(len(seq2)+1)]

    for i in range(1,len(seq1)+1):
        curr = [i*gap_penalty] + [0]*len(seq2)

        for j in range(1,len(seq2)+1):
            match = prev[j-1] + (match_score if seq1[i-1]==seq2[j-1] else mismatch_penalty)
            delete = prev[j] + gap_penalty
            insert = curr[j-1] + gap_penalty

            curr[j] = max(match,delete,insert)

        prev = curr

    return prev

def hirschberg(seq1, seq2):
    if len(seq1) == 0:
        return "-"*len(seq2), seq2
    if len(seq2) == 0:
        return seq1, "-"*len(seq1)

    if len(seq1) == 1 or len(seq2) == 1:
        return needleman_wunsch(seq1, seq2)

    mid = len(seq1)//2

    scoreL = nw_score(seq1[:mid], seq2)
    scoreR = nw_score(seq1[mid:][::-1], seq2[::-1])

    m = len(seq2)
    split = max(range(m+1), key=lambda j: scoreL[j] + scoreR[m-j])

    left1, left2 = hirschberg(seq1[:mid], seq2[:split])
    right1, right2 = hirschberg(seq1[mid:], seq2[split:])

    return left1 + right1, left2 + right2

def needleman_wunsch(seq1, seq2):
    match, mismatch, gap = 1, -1, -1
    n, m = len(seq1), len(seq2)

    dp = [[0]*(m+1) for _ in range(n+1)]

    for i in range(n+1):
        dp[i][0] = i * gap
    for j in range(m+1):
        dp[0][j] = j * gap

    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = dp[i-1][j] + gap
            insert = dp[i][j-1] + gap
            dp[i][j] = max(match_score, delete, insert)

    i, j = n, m
    a1, a2 = "", ""

    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + (
            match if seq1[i-1] == seq2[j-1] else mismatch
        ):
            a1 += seq1[i-1]
            a2 += seq2[j-1]
            i -= 1
            j -= 1
        elif i > 0:
            a1 += seq1[i-1]
            a2 += "-"
            i -= 1
        else:
            a1 += "-"
            a2 += seq2[j-1]
            j -= 1

    return a1[::-1], a2[::-1]


def global_alignment(seq1, seq2):
    align1, align2 = hirschberg(seq1, seq2)
    score = sum(
        1 if a == b else -1 if a != "-" and b != "-" else -1
        for a, b in zip(align1, align2)
    )
    return score, align1 + "\n" + align2



print(global_alignment("atgat","atgt"))


