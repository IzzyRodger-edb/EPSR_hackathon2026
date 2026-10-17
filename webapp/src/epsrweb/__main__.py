"""`python -m epsrweb` -- run the development server.

Jobs are started with multiprocessing's "spawn" method, so the reloader is off
by default: a reload would orphan a running refinement.
"""

from __future__ import annotations

import argparse
import logging
import os


def main() -> None:
    parser = argparse.ArgumentParser(prog="epsrweb", description=__doc__)
    parser.add_argument("--host", default=os.environ.get("EPSRWEB_HOST", "127.0.0.1"))
    parser.add_argument("--port", type=int, default=int(os.environ.get("EPSRWEB_PORT", "8000")))
    parser.add_argument("--log-level", default="info")
    args = parser.parse_args()

    logging.basicConfig(
        level=args.log_level.upper(),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    import uvicorn

    from .app import create_app

    uvicorn.run(create_app(), host=args.host, port=args.port, log_level=args.log_level)


if __name__ == "__main__":
    main()
