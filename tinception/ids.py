"""Subject identifiers.

Every subject gets one BIDS label ``sub-<sitelabel><sanitised original id>``. The site prefix is
mandatory because sites reuse IDs (e.g. 32 IDs appear in both UIUC and WHASC for different people).
"""

from __future__ import annotations

import re

from . import config

_NON_ALNUM = re.compile(r"[^A-Za-z0-9]")


def sanitise(original_id: str) -> str:
    """Strip everything BIDS disallows in a label (only [A-Za-z0-9] is allowed)."""
    return _NON_ALNUM.sub("", str(original_id).removeprefix("sub-"))


def bids_id(site: str, original_id: str) -> str:
    """Return ``sub-<sitelabel><id>`` for a site key from config/sites.yaml."""
    label = config.sites()[site]["label"]
    clean = sanitise(original_id)
    if not clean:
        raise ValueError(f"empty ID after sanitising {original_id!r} ({site})")
    return f"sub-{label}{clean}"


def canonical_original_id(site: str, original_id) -> str:
    """Normalise known per-site quirks so phenotype IDs match image/FreeSurfer IDs."""
    sid = str(original_id).strip()
    if site == "tinmeg" and sid.isdigit():
        sid = sid.zfill(4)  # images/FreeSurfer use 0539, spreadsheet has 539
    return sid
