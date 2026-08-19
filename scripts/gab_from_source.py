#!/usr/bin/env python3
"""Run the checked-out builder source inside a workflow-managed environment."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

main = importlib.import_module("genome_annotation_builder.cli").main


if __name__ == "__main__":
    raise SystemExit(main())
