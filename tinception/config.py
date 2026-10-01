"""Configuration loading: paths per environment, analysis parameters and the site registry."""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = REPO_ROOT / "config"


def _load(name: str) -> dict:
    with open(CONFIG_DIR / name) as f:
        return yaml.safe_load(f)


@dataclass(frozen=True)
class Paths:
    root: Path
    rawdata: Path
    derivatives: Path
    freesurfer: Path
    phenotypes_raw: Path
    legacy_masters: Path
    results: Path
    sources: dict[str, Path]


def _resolve(p: str) -> Path:
    path = Path(p)
    return path if path.is_absolute() else REPO_ROOT / path


@lru_cache
def paths(env: str | None = None) -> Paths:
    env = env or os.environ.get("TINCEPTION_ENV", "local")
    cfg = _load("paths.yaml")[env]
    return Paths(
        root=_resolve(cfg["root"]),
        rawdata=_resolve(cfg["rawdata"]),
        derivatives=_resolve(cfg["derivatives"]),
        freesurfer=_resolve(cfg["freesurfer"]),
        phenotypes_raw=_resolve(cfg["phenotypes_raw"]),
        legacy_masters=_resolve(cfg["legacy_masters"]),
        results=_resolve(cfg["results"]),
        sources={k: _resolve(v) for k, v in (cfg.get("sources") or {}).items()},
    )


@lru_cache
def analysis() -> dict:
    return _load("analysis.yaml")


@lru_cache
def sites() -> dict[str, dict]:
    return _load("sites.yaml")["sites"]
