"""
Syntactic validation of H5MD ``unit`` attribute strings.

S4.2 of the specification requires a ``unit`` attribute on every dataset holding
a physical quantity, and recommends the notation of the H5MD ``units`` module.
This module checks that notation:

    unit_string := [ prefactor SPACE ] factor { SPACE factor }
    factor      := symbol [ exponent ]
    exponent    := ( "+" | "-" ) digit+
    prefactor   := decimal [ exponent ]      e.g. "10+5", "4.184"

Exponents are always written with an explicit sign, which is what distinguishes
``nm+3`` (a volume) from a symbol that merely ends in a digit.

Used by:
    Validation 3: internal consistency
"""

from __future__ import annotations

import re

_EXPONENT = r"(?:[+-]\d+)"
_PREFACTOR_RE = re.compile(rf"^\d+(?:\.\d+)?{_EXPONENT}?$")
_FACTOR_RE = re.compile(rf"^[A-Za-z_%°µΩ]+{_EXPONENT}?$")

#: Symbols used by the H5MD units module plus the non-SI units the
#: specification explicitly admits (S8 note: degree is accepted for use with SI).
KNOWN_SYMBOLS = {
    # SI base
    "m",
    "kg",
    "s",
    "A",
    "K",
    "mol",
    "cd",
    # SI derived / accepted
    "Hz",
    "N",
    "Pa",
    "J",
    "W",
    "C",
    "V",
    "F",
    "Ohm",
    "S",
    "Wb",
    "T",
    "H",
    "lm",
    "lx",
    "Bq",
    "Gy",
    "Sv",
    "kat",
    "rad",
    "sr",
    "degree",
    "deg",
    "bar",
    "u",
    "e",
    "eV",
    "L",
    "l",
    "cal",
    "amu",
    "Da",
    # common prefixed forms appearing in MD output
    "nm",
    "pm",
    "fm",
    "um",
    "mm",
    "cm",
    "km",
    "ps",
    "fs",
    "ns",
    "us",
    "ms",
    "kJ",
    "kcal",
    "kPa",
    "MPa",
    "GPa",
    "kN",
    "mN",
    "nN",
    "pN",
    "g",
    "mg",
    "ng",
}


def parse(unit: str):
    """Split a unit string into (prefactor_token, [factor_tokens]).

    Returns ``None`` when the string is syntactically invalid.
    """
    if unit is None:
        return None
    tokens = unit.strip().split()
    if not tokens:
        # The empty string is legal: it marks a dimensionless quantity.
        return ("", [])
    prefactor = ""
    if _PREFACTOR_RE.match(tokens[0]) and not _FACTOR_RE.match(tokens[0]):
        prefactor, tokens = tokens[0], tokens[1:]
    elif _PREFACTOR_RE.match(tokens[0]):
        # Purely numeric leading token, e.g. "4.184" in "4.184 kJ mol-1".
        prefactor, tokens = tokens[0], tokens[1:]
    for tok in tokens:
        if not _FACTOR_RE.match(tok):
            return None
    return (prefactor, tokens)


def is_valid(unit: str) -> bool:
    """True when ``unit`` is syntactically well formed."""
    return parse(unit) is not None


def symbol_of(factor: str) -> str:
    """Strip the exponent from a factor token: ``nm-2`` -> ``nm``."""
    return re.sub(rf"{_EXPONENT}$", "", factor)


def unknown_symbols(unit: str):
    """Symbols in ``unit`` that are not in :data:`KNOWN_SYMBOLS`."""
    parsed = parse(unit)
    if parsed is None:
        return []
    _, factors = parsed
    return [s for s in (symbol_of(f) for f in factors) if s not in KNOWN_SYMBOLS]
