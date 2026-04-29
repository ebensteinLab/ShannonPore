"""Parse Region-Of-Interest definitions from BED-like files or text areas."""

from __future__ import annotations

import logging
import re

import pandas as pd

logger = logging.getLogger(__name__)


def parse_roi_file_flexible(path: str) -> pd.DataFrame | None:
    """Parse a BED-like ROI file. Accepts header with chr/start/end (and
    optional target) OR no header (col1=chr col2=start col3=end col4=target).

    Returns None on parse failure. Returns DataFrame with columns
    [chr, start, end, target] otherwise.
    """
    try:
        df_raw = pd.read_csv(path, sep="\t", comment="#", header=None)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        logger.warning("parse_roi_file_flexible(%s): %s", path, exc)
        return None

    if df_raw.shape[1] < 3:
        return None

    first_row = df_raw.iloc[0].astype(str).str.lower()
    header_like = (
        first_row.iloc[0] in ["chr", "chrom", "chromosome"]
        and first_row.iloc[1] in ["start", "pos", "position"]
        and first_row.iloc[2] in ["end", "stop"]
    )

    try:
        if header_like:
            cols = list(df_raw.iloc[0])
            df = df_raw.iloc[1:].copy()
            df.columns = cols

            rename: dict[str, str] = {}
            for c in df.columns:
                lc = str(c).lower()
                if lc in ["chr", "chrom", "chromosome"]:
                    rename[c] = "chr"
                elif lc in ["start", "pos", "position"]:
                    rename[c] = "start"
                elif lc in ["end", "stop"]:
                    rename[c] = "end"
                elif lc in ["target", "name", "gene", "id"]:
                    rename[c] = "target"
            df = df.rename(columns=rename)

            if "target" not in df.columns:
                df["target"] = (
                    df["chr"].astype(str)
                    + ":" + df["start"].astype(str)
                    + "-" + df["end"].astype(str)
                )

            df = df[["chr", "start", "end", "target"]]
        else:
            df = df_raw.iloc[:, :4].copy()
            df.columns = ["chr", "start", "end", "target"][: df.shape[1]]
            if "target" not in df.columns:
                df["target"] = (
                    df["chr"].astype(str)
                    + ":" + df["start"].astype(str)
                    + "-" + df["end"].astype(str)
                )
            df = df[["chr", "start", "end", "target"]]

        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        return df
    except (ValueError, KeyError) as exc:
        logger.warning("parse_roi_file_flexible(%s) malformed: %s", path, exc)
        return None


def parse_roi_text_area(text_value: str) -> pd.DataFrame | None:
    """Parse pasted ROI in BED-like format: `chr\tstart\tend\t[name]`."""
    if not text_value or text_value.strip() == "":
        return None

    lines = [
        ln.strip()
        for ln in text_value.strip().splitlines()
        if ln.strip() and not ln.strip().startswith("#")
    ]
    if not lines:
        return None

    rows = []
    for ln in lines:
        parts = re.split(r"\s+", ln)
        if len(parts) < 3:
            return None
        chrom = parts[0]
        try:
            start = int(parts[1])
            end = int(parts[2])
        except ValueError:
            return None
        if end <= start:
            return None
        target = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
        rows.append({"chr": chrom, "start": start, "end": end, "target": target})

    return pd.DataFrame(rows)
