"""Mellor: an LLM agent for EPSR refinements.

The main entry points in one place::

    from mellor import Mellor, SystemSpec, make_backend, run_spec

``SystemSpec`` is the specification the model, the library or a user produces;
``run_spec`` executes one with the EPSR engine; ``make_backend`` picks the language
model (torch and transformers are imported only when such a backend is used);
``Mellor`` is the agent. See docs/architecture.md and docs/integration.md.
"""

from ._version import __version__
from .llm.agent import Mellor
from .llm.backend import make_backend
from .runner import run_spec
from .spec import SystemSpec

__all__ = ["__version__", "Mellor", "SystemSpec", "make_backend", "run_spec"]
