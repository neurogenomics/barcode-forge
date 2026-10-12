"""Reverse-complement palindromes.

Palindromic stems self-anneal (primer dimers, hairpins) and are typical
restriction sites (EcoRI GAATTC, HindIII AAGCTT).
"""

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule

_PAIR = {"A": "T", "T": "A", "C": "G", "G": "C"}


def longest_palindrome(seq: str) -> int:
    """Length of the longest substring equal to its reverse complement (always even)."""
    best = 0
    for centre in range(1, len(seq)):
        left, right = centre - 1, centre
        while left >= 0 and right < len(seq) and _PAIR[seq[left]] == seq[right]:
            left, right = left - 1, right + 1
        best = max(best, right - left - 1)
    return best


def check(seq: str, cfg: Config) -> bool:
    return cfg.max_palindrome is None or longest_palindrome(seq) <= cfg.max_palindrome


RULE = Rule(name="palindrome", cost=2, check=check)
