"""DNA alphabet and basic sequence helpers."""

DNA = "ACGT"

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def validate(seq: str) -> str:
    """Return ``seq`` uppercased, or raise ValueError if it is not non-empty ACGT."""
    upper = seq.upper()
    if not upper or set(upper) - set(DNA):
        raise ValueError(f"not a DNA sequence (ACGT only): {seq!r}")
    return upper


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    if not seq:
        raise ValueError("GC fraction of an empty sequence is undefined")
    return (seq.count("G") + seq.count("C")) / len(seq)
