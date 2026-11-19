"""Read simulation for barcode layouts."""

from .error_model import (
    DEFAULT_DRATE,
    DEFAULT_IRATE,
    DEFAULT_SRATE,
    reverse_complement,
    simulate_phred,
    simulate_read,
)

__all__ = [
    "DEFAULT_SRATE",
    "DEFAULT_IRATE",
    "DEFAULT_DRATE",
    "reverse_complement",
    "simulate_read",
    "simulate_phred",
]
