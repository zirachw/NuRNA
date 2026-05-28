from typing import Generator


def can_pair(a: str, b: str) -> bool:
    return (a, b) in {("A","U"),("U","A"),("G","C"),("C","G"),("G","U"),("U","G")}


def build_dp_matrix(seq: str) -> list[list[int]]:
    n = len(seq)
    dp = [[0] * n for _ in range(n)]
    for length in range(1, n):
        for i in range(n - length):
            j = i + length
            # gamma[i][j] = max(unpaired_i, unpaired_j, paired(i,j) if compatible, max_k bifurcation)
            dp[i][j] = max(
                dp[i + 1][j],
                dp[i][j - 1],
                (dp[i + 1][j - 1] + 1) if can_pair(seq[i], seq[j]) and j - i > 3 else 0,
                max((dp[i][k] + dp[k + 1][j]) for k in range(i, j)),
            )
    return dp


def build_dp_matrix_steps(seq: str) -> Generator[list[list[int]], None, None]:
    n = len(seq)
    dp = [[0] * n for _ in range(n)]
    for length in range(1, n):
        for i in range(n - length):
            j = i + length
            # gamma[i][j] = max(unpaired_i, unpaired_j, paired(i,j) if compatible, max_k bifurcation)
            dp[i][j] = max(
                dp[i + 1][j],
                dp[i][j - 1],
                (dp[i + 1][j - 1] + 1) if can_pair(seq[i], seq[j]) and j - i > 3 else 0,
                max((dp[i][k] + dp[k + 1][j]) for k in range(i, j)),
            )
            yield [row[:] for row in dp]


def traceback(
    dp: list[list[int]],
    seq: str,
    i: int,
    j: int,
    pairs: set[tuple[int, int]],
) -> None:
    if i >= j:
        return
    # i is unpaired: optimal score unchanged when i is skipped
    if dp[i][j] == dp[i + 1][j]:
        traceback(dp, seq, i + 1, j, pairs)
    # j is unpaired: optimal score unchanged when j is skipped
    elif dp[i][j] == dp[i][j - 1]:
        traceback(dp, seq, i, j - 1, pairs)
    # i and j are paired: score came from gamma[i+1][j-1] + 1
    elif can_pair(seq[i], seq[j]) and j - i > 3 and dp[i][j] == dp[i + 1][j - 1] + 1:
        pairs.add((i, j))
        traceback(dp, seq, i + 1, j - 1, pairs)
    # bifurcation: find the split point k where gamma[i][k] + gamma[k+1][j] = gamma[i][j]
    else:
        for k in range(i + 1, j):
            if dp[i][j] == dp[i][k] + dp[k + 1][j]:
                traceback(dp, seq, i, k, pairs)
                traceback(dp, seq, k + 1, j, pairs)
                break


def predict_structure(seq: str) -> str:
    seq = seq.upper()
    dp = build_dp_matrix(seq)
    pairs: set[tuple[int, int]] = set()
    traceback(dp, seq, 0, len(seq) - 1, pairs)
    result = ["."] * len(seq)
    for i, j in pairs:
        result[i] = "("
        result[j] = ")"
    return "".join(result)
