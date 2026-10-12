"""Sequence-Levenshtein distance for inline barcodes.

Buschmann T, Bystrykh LV (2013). Levenshtein error-correcting barcodes for
multiplexed DNA sequencing. BMC Bioinformatics 14:272. doi:10.1186/1471-2105-14-272

An inline barcode is followed by downstream sequence, so an indel shifts the read:
a deletion pulls a downstream base into the barcode window and an insertion pushes
the last barcode base out of it. Truncating or extending either sequence at its
3' end is therefore free. This is the Levenshtein DP matrix read out as the
minimum over its last row and last column instead of its corner.

It is not a metric: a sequence is 0 away from any extension of itself, and the
triangle inequality can fail.
"""


def distance(a: str, b: str, *, score_cutoff: int | None = None) -> int:
    """Sequence-Levenshtein distance; returns ``score_cutoff + 1`` if it exceeds the cutoff."""
    cap = len(a) + len(b) + 1 if score_cutoff is None else score_cutoff + 1
    prev = list(range(len(b) + 1))
    best = prev[-1]  # minimum of the last column so far
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
        best = min(best, cur[-1])
        # Row minima never decrease, so no later cell can beat min(best, cap).
        if min(cur) >= min(best, cap):
            return min(best, cap)
    return min(best, min(prev), cap)
