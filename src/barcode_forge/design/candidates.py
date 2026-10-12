"""Seeded stream of unique random candidate sequences."""

from collections.abc import Callable
from random import Random

from barcode_forge.sequence import DNA

# At or below this space size every sequence is enumerated, so a small design
# sees the whole space (4**10 = 2**20 sequences).
EXHAUSTIVE_MAX = 2**20
# Random mode gives up after this many draws per requested sequence, so a strict
# filter can never hang the generator.
DRAWS_PER_LIMIT = 50


def _decode(index: int, length: int) -> str:
    chars = []
    for _ in range(length):
        index, base = divmod(index, 4)
        chars.append(DNA[base])
    return "".join(chars)


def generate(
    length: int,
    limit: int,
    rng: Random,
    accept: Callable[[str], bool] = lambda s: True,
) -> list[str]:
    """Return up to ``limit`` unique sequences of ``length`` passing ``accept``, in random order.

    Fewer than ``limit`` are returned only if the space (or the draw budget in
    random mode) runs out; callers decide whether that is fatal.
    """
    if length < 1:
        raise ValueError("length must be >= 1")
    if limit < 0:
        raise ValueError("limit must be >= 0")
    out: list[str] = []
    if limit == 0:
        return out
    space = 4**length
    if space <= EXHAUSTIVE_MAX:
        order = list(range(space))
        rng.shuffle(order)
        for index in order:
            seq = _decode(index, length)
            if accept(seq):
                out.append(seq)
                if len(out) == limit:
                    break
        return out
    seen: set[str] = set()
    for _ in range(DRAWS_PER_LIMIT * limit):
        seq = "".join(rng.choices(DNA, k=length))
        if seq in seen:
            continue
        seen.add(seq)
        if accept(seq):
            out.append(seq)
            if len(out) == limit:
                break
    return out
