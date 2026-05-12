from __future__ import annotations

from pathlib import Path
from typing import Literal

import polars as pl
import polars.selectors as cs

from .parameters import (
    _DIA_FP_CONSTANT_ION_COLUMNS,
    _DIA_RENAME_FP_ION,
    _DIA_RENAME_DIANN_ION,
    _DIA_RENAME_DIANN_PROTEIN,
)

type Annotation = dict[str, dict[str, str]]
type DataType = Literal["ion", "trp"]


def _read_ion(path: Path, method: str) -> pl.DataFrame:
    match method:
        case "dda":
            dda_ion_df = pl.read_csv(path.joinpath("combined_ion.tsv"), separator="\t")
            dda_ion_df = _zero_fill_intensities(dda_ion_df)

            return dda_ion_df

        case "dia":
            annot = _read_experiment_annotation(path)

            dia_ion_df = pl.read_csv(
                path.joinpath("dia-quant-output/report.pr_matrix.tsv"),
                separator="\t",
            )
            fp_ion_df = pl.read_csv(path.joinpath("ion.tsv"), separator="\t").select(
                _DIA_FP_CONSTANT_ION_COLUMNS
            )

            dia_ion_df = _rename_dia_columns(dia_ion_df, annot)
            dia_ion_df = _join_dia_ion_metadata(dia_ion_df, fp_ion_df)
            dia_ion_df = _zero_fill_intensities(dia_ion_df)

            return dia_ion_df

        case _:
            raise ValueError("Input error.")


def _read_trp(path: Path, method: str) -> pl.DataFrame:
    match method:
        case "dda":
            dda_trp_df = pl.read_csv(path.joinpath("combined_protein.tsv"), separator="\t")
            dda_trp_df = _zero_fill_intensities(dda_trp_df)

            return dda_trp_df

        case "dia":
            annot = _read_experiment_annotation(path)

            dia_trp_df = pl.read_csv(
                path.joinpath("dia-quant-output/report.pg_matrix.tsv"),
                separator="\t",
            )

            dia_trp_df = _rename_dia_columns(dia_trp_df, annot, "trp")
            dia_trp_df = _zero_fill_intensities(dia_trp_df)

            return dia_trp_df

        case _:
            raise ValueError("Input error.")


def _zero_fill_intensities(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(cs.ends_with("Intensity").fill_null(0))


def _read_experiment_annotation(path: Path) -> Annotation:
    df = pl.read_csv(
        path.joinpath("experiment_annotation.tsv"),
        separator="\t",
    )

    annot: Annotation = {}
    for row in df.iter_rows():
        file, sample, sample_name, condition, replicate = row
        annot[str(file).rsplit(".", maxsplit=1)[0]] = {
            "Sample": str(sample),
            "Sample Name": str(sample_name),
            "Condition": str(condition),
            "Replicate": str(replicate),
        }

    return annot


def _file_to_sample_name(annot: Annotation) -> dict[str, str]:
    return {file: data.get("Sample Name", "") for file, data in annot.items()}


def _rename_dia_columns(
    df: pl.DataFrame,
    annot: Annotation,
    data_type: DataType = "ion",
) -> pl.DataFrame:
    other_cols = _DIA_RENAME_DIANN_ION
    if data_type != "ion":
        other_cols = _DIA_RENAME_DIANN_PROTEIN

    rename: dict[str, str] = {}
    for file, sample in _file_to_sample_name(annot).items():
        for col in df.columns:
            if file in col:
                rename[col] = f"{sample} Intensity"

    return df.rename(rename).rename(other_cols)


def _join_dia_ion_metadata(dia_df: pl.DataFrame, ion_df: pl.DataFrame) -> pl.DataFrame:
    dia_df = dia_df.with_columns(
        pl.concat_str("Protein ID", "Peptide Sequence").alias("Unique ID")
    )

    ion_df = (
        ion_df.with_columns(
            pl.concat_str("Protein ID", "Peptide Sequence").alias("Unique ID")
        )
        .drop("Protein ID", "Peptide Sequence")
        .unique("Unique ID")
    )

    ion_df = ion_df.rename(_DIA_RENAME_FP_ION)

    return dia_df.join(ion_df, on="Unique ID", how="left")
