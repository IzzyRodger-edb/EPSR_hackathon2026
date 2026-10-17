"""Exceptions raised by the EPSR engine.

Every engine error is an :class:`EpsrError`. Each subclass also derives from the
matching built-in exception, so ``except OSError`` or ``except ValueError`` work too.

The Rust extension (``epsr._core``) raises these classes directly.
"""

from __future__ import annotations

__all__ = ["EpsrError", "EpsrIOError", "EpsrInputError", "EpsrUnsupportedError"]


class EpsrError(Exception):
    """Base class of all errors reported by the EPSR engine."""


class EpsrIOError(EpsrError, OSError):
    """A file could not be opened, read or written (Rust ``EpsrError::Io``)."""


class EpsrInputError(EpsrError, ValueError):
    """An input file or override line is malformed (Rust ``EpsrError::Input``)."""


class EpsrUnsupportedError(EpsrError, NotImplementedError):
    """The input needs a feature on which the original Fortran itself fails.

    Rust ``EpsrError::Unsupported``. The list of such inputs is in
    ``rust/docs/FORTRAN_MAP.md`` ("Still unsupported").
    """
