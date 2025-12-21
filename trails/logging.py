# logging_utils.py
from __future__ import annotations

import logging
import logging.handlers
import os
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from typing import Optional

# Context variables to correlate logs across modules
_run_id: ContextVar[str] = ContextVar("run_id", default="-")
_case: ContextVar[str] = ContextVar("case", default="-")
_year: ContextVar[str] = ContextVar("year", default="-")
_depth: ContextVar[str] = ContextVar("depth", default="-")


class TrailsContextFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        record.run_id = _run_id.get()
        record.case = _case.get()
        record.year = _year.get()
        record.depth = _depth.get()
        return True


def configure_trails_logging(
    file_level: int = logging.DEBUG,
    filename: str = "trails.log",
    also_console: bool = False,
    console_level: int = logging.WARNING,
    mode: str = "w",
    debug: bool = False,
) -> Path:
    """
    Configure root logging so all modules using logging.getLogger(__name__)
    write to the same file in the current working directory.

    - File handler always writes to `filename`.
    - Console handler (optional) is kept at WARNING by default (notebook-friendly).
    - Idempotent across notebook reruns: updates or replaces existing Trails handlers.
    """
    log_path = Path(os.getcwd()) / filename
    root = logging.getLogger()
    root.setLevel(min(file_level, console_level) if also_console else file_level)

    fmt = (
        "%(asctime)s | %(levelname)-7s | %(name)s | "
        "run=%(run_id)s case=%(case)s year=%(year)s depth=%(depth)s | "
        "%(message)s"
    )
    formatter = logging.Formatter(fmt)

    # Helper to detect "our" handlers
    def _is_trails_file_handler(h: logging.Handler) -> bool:
        return isinstance(h, logging.FileHandler) and getattr(
            h, "baseFilename", ""
        ) == str(log_path)

    def _is_console_handler(h: logging.Handler) -> bool:
        # StreamHandler that is not a FileHandler
        return isinstance(h, logging.StreamHandler) and not isinstance(
            h, logging.FileHandler
        )

    # --- Ensure exactly one file handler for this log_path ---
    file_handlers = [h for h in root.handlers if _is_trails_file_handler(h)]
    if file_handlers:
        fh = file_handlers[0]
        fh.setLevel(file_level)
        fh.setFormatter(formatter)
        # Ensure filter present (avoid duplicates)
        fh.filters = [f for f in fh.filters if not isinstance(f, TrailsContextFilter)]
        fh.addFilter(TrailsContextFilter())
        # Remove extra duplicates if any
        for extra in file_handlers[1:]:
            root.removeHandler(extra)
    else:
        fh = logging.FileHandler(log_path, mode=mode, encoding="utf-8")
        fh.setLevel(file_level)
        fh.setFormatter(formatter)
        fh.addFilter(TrailsContextFilter())
        root.addHandler(fh)

    # --- Console handler control (notebook) ---
    if also_console:
        # Prefer to *reuse* one existing console handler; otherwise create one
        console_handlers = [h for h in root.handlers if _is_console_handler(h)]
        if console_handlers:
            ch = console_handlers[0]
            ch.setLevel(console_level)
            ch.setFormatter(formatter)
            ch.filters = [
                f for f in ch.filters if not isinstance(f, TrailsContextFilter)
            ]
            ch.addFilter(TrailsContextFilter())
            # Remove extra duplicates if any
            for extra in console_handlers[1:]:
                root.removeHandler(extra)
        else:
            ch = logging.StreamHandler()
            ch.setLevel(console_level)
            ch.setFormatter(formatter)
            ch.addFilter(TrailsContextFilter())
            root.addHandler(ch)
    else:
        # Remove console handlers if user does not want console output
        for h in list(root.handlers):
            if _is_console_handler(h):
                root.removeHandler(h)

    # Keep noisy libs from drowning your file unless you explicitly want them
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)

    if debug:
        root.info(
            "Logging configured. File=%s (level=%s) Console=%s",
            str(log_path),
            logging.getLevelName(file_level),
            logging.getLevelName(console_level) if also_console else "OFF",
        )

    return log_path


@contextmanager
def trails_log_context(
    run_id: Optional[str] = None,
    case: Optional[str] = None,
    year: Optional[int] = None,
    depth: Optional[int] = None,
):
    tokens = []
    try:
        if run_id is not None:
            tokens.append((_run_id, _run_id.set(str(run_id))))
        if case is not None:
            tokens.append((_case, _case.set(str(case))))
        if year is not None:
            tokens.append((_year, _year.set(str(year))))
        if depth is not None:
            tokens.append((_depth, _depth.set(str(depth))))
        yield
    finally:
        for var, tok in reversed(tokens):
            var.reset(tok)
