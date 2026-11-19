"""Article-style incremental noisy channel for barcode read simulation.

Per reference position, draw one of {match, substitution, deletion, insertion}
with unequal probabilities. Insertion emits a base and does **not** advance the
reference cursor. No pad/truncate — read length is whatever the walk emits.
"""

from __future__ import annotations

import random
from typing import Tuple

BASES = ("A", "C", "G", "T")
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

DEFAULT_SRATE = 0.04
DEFAULT_IRATE = 0.05
DEFAULT_DRATE = 0.14


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def simulate_read(
    barcode: str,
    srate: float = DEFAULT_SRATE,
    irate: float = DEFAULT_IRATE,
    drate: float = DEFAULT_DRATE,
    rng: random.Random | None = None,
) -> Tuple[str, int, int, int]:
    """Walk ``barcode`` and return ``(read, n_sub, n_ins, n_del)``.

    At each step while the reference cursor has not finished:
      match         with prob 1 - srate - irate - drate  → emit ref base, advance
      substitution  with prob srate                     → emit other base, advance
      deletion      with prob drate                     → emit nothing, advance
      insertion     with prob irate                     → emit random base, stay
    """
    if rng is None:
        rng = random.Random()
    total = srate + irate + drate
    if total < 0 or total > 1:
        raise ValueError(f"srate+irate+drate must be in [0, 1], got {total}")

    p_match = 1.0 - total
    seq = barcode.upper()
    out: list[str] = []
    n_sub = n_ins = n_del = 0
    i = 0
    n = len(seq)

    while i < n:
        u = rng.random()
        if u < p_match:
            out.append(seq[i])
            i += 1
        elif u < p_match + srate:
            choices = [b for b in BASES if b != seq[i]]
            out.append(rng.choice(choices))
            n_sub += 1
            i += 1
        elif u < p_match + srate + drate:
            n_del += 1
            i += 1
        else:
            out.append(rng.choice(BASES))
            n_ins += 1

    return "".join(out), n_sub, n_ins, n_del


def simulate_phred(
    length: int,
    rng: random.Random,
    mean_q: float = 30.0,
    sd_q: float = 3.0,
    qmin: int = 10,
    qmax: int = 40,
) -> str:
    """Phred+33 quality string; i.i.d. Normal(mean_q, sd_q) clipped to [qmin, qmax]."""
    if length < 0:
        raise ValueError("length must be ≥ 0")
    chars = []
    for _ in range(length):
        q = int(round(rng.gauss(mean_q, sd_q)))
        q = max(qmin, min(qmax, q))
        chars.append(chr(q + 33))
    return "".join(chars)
