from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Final, Literal

import polars as pl
import polars.selectors as cs
import scipy as sp

from .functions import _add_log2_column, _add_neg_log10_column
from .parameters import _FLIPPR_COMBINE_KEY

type CombineLevel = Literal["CUT SITE", "PEPTIDE", "MODIFIED PEPTIDE"]

COMB_NAME_COLUMN: Final[dict[CombineLevel, str]] = {
    "CUT SITE": "Cut Site ID",
    "PEPTIDE": "Peptide Sequence",
    "MODIFIED PEPTIDE": "Modified Sequence",
}


def _combine_pvalues(p_values: Sequence[float]) -> float:
    return float(sp.stats.combine_pvalues(p_values)[1])


def combine_by(df: pl.DataFrame, by: CombineLevel, fc: str) -> pl.DataFrame:
    combined = (
        df.group_by(["Protein ID", COMB_NAME_COLUMN[by]], maintain_order=True)
        .agg(
            pl.col(_FLIPPR_COMBINE_KEY[by]).first(),
            cs.by_name("P-value", "Adj. P-value", "CV", fc).filter(
                pl.col("T-test").sign() == pl.col("T-test").sign().sum().sign()
            ),
        )
        .with_columns(
            pl.when(pl.col("P-value").list.len().gt(0))
            .then(pl.col("P-value"))
            .otherwise([1.0])
            .alias("P-value"),
            pl.when(pl.col("Adj. P-value").list.len().gt(0))
            .then(pl.col("Adj. P-value"))
            .otherwise([1.0])
            .alias("Adj. P-value"),
            pl.when(pl.col("CV").list.len().gt(0))
            .then(pl.col("CV"))
            .otherwise([0.0])
            .alias("CV"),
            pl.when(pl.col(fc).list.len().gt(0))
            .then(pl.col(fc))
            .otherwise([0.0])
            .alias(fc),
        )
        .with_columns(
            pl.col("P-value")
            .map_elements(_combine_pvalues, return_dtype=pl.Float64)
            .alias("P-value"),
            pl.col("Adj. P-value")
            .map_elements(_combine_pvalues, return_dtype=pl.Float64)
            .alias("Adj. P-value"),
            pl.col("CV").list.max()
            .alias("CV"),
            pl.col(fc).list.median()
            .alias(fc),
        )
    )

    combined = _add_log2_column(combined, fc)
    combined = _add_neg_log10_column(combined, "P-value")
    combined = _add_neg_log10_column(combined, "Adj. P-value")

    return combined


def summary_by(df: pl.DataFrame, by: str, fc: str, rcParams: Mapping[str, Any]) -> pl.DataFrame:
    prot_fc_sig = rcParams.get("protein.fc_sig_thresh", 1.0)
    prot_pv_sig = rcParams.get("protein.pval_sig_thresh", 0.01)
    prot_apv_sig = rcParams.get("protein.adj_pval_sig_thresh", 0.05)

    val = (
        df.group_by("Protein ID", maintain_order=True)
        .agg(pl.col("CV").filter(pl.col("-Log10 P-value").gt(0)).len())
        .rename({"CV": f"No. of Valid {by}"})
    )

    sig = (
        df.group_by("Protein ID", maintain_order=True)
        .agg(
            pl.col("CV").filter(
                (pl.col(f"Log2 {fc}").abs().ge(prot_fc_sig))
                & (pl.col("P-value").le(prot_pv_sig))
            ).len()
        )
        .rename({"CV": f"No. of Significant {by} (P-value)"})
    )

    sigsig = (
        df.group_by("Protein ID", maintain_order=True)
        .agg(
            pl.col("CV").filter(
                (pl.col(f"Log2 {fc}").abs().ge(prot_fc_sig))
                & (pl.col("Adj. P-value").le(prot_apv_sig))
            ).len()
        )
        .rename({"CV": f"No. of Significant {by} (Adj. P-value)"})
    )

    return val.join(sig, on="Protein ID").join(sigsig, on="Protein ID")
