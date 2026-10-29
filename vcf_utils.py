"""Shared helpers for reading single-nucleotide variants from VCF files.

Both ``snv_substitution_spectrum.py`` and ``variant_concordance_analysis.py``
build on the same two primitives: pulling single-base substitutions out of a
(gzipped) VCF file and collapsing each substitution onto one of the six
strand-symmetric SNV classes. Keeping that logic in one place avoids drift
between the two analyses.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import NamedTuple

# Strand-symmetric collapse of the 12 possible single-base substitutions
# onto the 6 canonical classes (pyrimidine-based notation), i.e. the same
# convention used by the COSMIC single-base-substitution (SBS) signatures.
_SNV_CLASSES = {
    "C>A": "C>A", "G>T": "C>A",
    "C>G": "C>G", "G>C": "C>G",
    "C>T": "C>T", "G>A": "C>T",
    "T>A": "T>A", "A>T": "T>A",
    "T>C": "T>C", "A>G": "T>C",
    "T>G": "T>G", "A>C": "T>G",
}

SNV_CLASS_ORDER = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


class Variant(NamedTuple):
    """A minimal, hashable representation of a single VCF record."""

    chrom: str
    pos: str
    ref: str
    alt: str


def classify_snv(ref: str, alt: str) -> str | None:
    """Collapse a single-base substitution onto one of the 6 SNV classes.

    Returns ``None`` if ``ref``/``alt`` is not a single-base substitution
    (e.g. an indel or multi-allelic record), so callers can filter those
    out with a simple truthiness check.
    """
    return _SNV_CLASSES.get(f"{ref}>{alt}")


def normalize_chromosome(chrom: str) -> str:
    """Strip a leading ``chr``/``Chr`` prefix and uppercase mitochondrial names.

    This makes chromosome names comparable across VCFs that were called
    with different naming conventions (e.g. ``chr1`` vs ``1``, ``chrM`` vs
    ``MT``).
    """
    normalized = chrom.replace("chr", "").replace("Chr", "")
    return "M" if normalized.upper() == "MT" else normalized


def _open_vcf(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_variants(vcf_path: Path, normalize_chrom: bool = False) -> list[Variant]:
    """Read every variant record from a (optionally gzipped) VCF file.

    Parameters
    ----------
    vcf_path:
        Path to a ``.vcf`` or ``.vcf.gz`` file.
    normalize_chrom:
        If True, chromosome names are passed through :func:`normalize_chromosome`
        before being stored, so that e.g. ``chr1`` and ``1`` are treated as
        the same chromosome when comparing variant sets across files.
    """
    variants: list[Variant] = []
    with _open_vcf(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _id, ref, alt = fields[:5]
            if normalize_chrom:
                chrom = normalize_chromosome(chrom)
            variants.append(Variant(chrom, pos, ref, alt))
    return variants


def read_snv_variants(vcf_path: Path, normalize_chrom: bool = False) -> list[Variant]:
    """Like :func:`read_variants`, but keeps only single-base substitutions."""
    return [
        v for v in read_variants(vcf_path, normalize_chrom=normalize_chrom)
        if len(v.ref) == 1 and len(v.alt) == 1
    ]
