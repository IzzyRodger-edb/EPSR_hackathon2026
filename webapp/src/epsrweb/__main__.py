"""`uv run epsrweb` (or `python -m epsrweb`) -- run the development server.

Jobs are started with multiprocessing's "spawn" method, so the reloader is off
by default: a reload would orphan a running refinement.

Ctrl+C stops the server in a second or two. The moment it is pressed, the
chat process is stopped (a question in flight fails at once) and the progress
streams close; the app's shutdown then stops the running refinements. Anything
still open after `SHUTDOWN_GRACE` seconds is cancelled.
"""

from __future__ import annotations

import argparse
import logging
import os
import threading

#: Seconds uvicorn waits for open requests after Ctrl+C before cancelling them.
#: A backstop: the progress streams and the chat end by themselves on shutdown.
SHUTDOWN_GRACE = 5


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

    app = create_app()

    class Server(uvicorn.Server):
        def handle_exit(self, sig, frame) -> None:
            super().handle_exit(sig, frame)
            # uvicorn waits for open requests before it runs the app's shutdown.
            # End the long ones now: the progress streams watch this flag, and
            # stopping the chat process makes a question in flight return at
            # once. In a thread, because this runs inside a signal handler.
            app.state.shutting_down.set()
            threading.Thread(target=app.state.jobs.chat.stop, name="epsrweb-chat-stop",
                             daemon=True).start()

    config = uvicorn.Config(app, host=args.host, port=args.port, log_level=args.log_level,
                            timeout_graceful_shutdown=SHUTDOWN_GRACE)
    try:
        Server(config).run()
    except KeyboardInterrupt:
        pass
    finally:
        # A second Ctrl+C makes uvicorn skip the app's shutdown. Stop the jobs here
        # too, or Python would wait for them before it exits.
        app.state.jobs.stop(terminate_running=True)


if __name__ == "__main__":
    main()
