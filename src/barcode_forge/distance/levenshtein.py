"""Levenshtein (edit) distance with unit costs for substitution, insertion and deletion.

``distance(a, b, *, score_cutoff=None)`` is rapidfuzz's C scorer itself (not a
wrapper) so ``design.pool`` can scan whole pools with ``rapidfuzz.process``.
"""

from rapidfuzz.distance.Levenshtein import distance

__all__ = ["distance"]
