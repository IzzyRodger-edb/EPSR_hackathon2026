"""A web application for the Mellor EPSR agent.

Upload measured structure-factor data, let the agent set up and run the Rust
EPSR engine, read the results as seaborn figures and a PDF report, and publish
chosen figures straight into the `figures/` folder of an Overleaf project.

The agent's side of the contract is `mellor/docs/integration.md`; the module
`mellor_adapter` is the only place that calls it.
"""

__version__ = "0.1.0"

__all__ = ["__version__", "create_app"]


def create_app(*args, **kwargs):
    """Lazy re-export, so `import epsrweb` stays cheap for tests and tooling."""
    from .app import create_app as factory

    return factory(*args, **kwargs)
