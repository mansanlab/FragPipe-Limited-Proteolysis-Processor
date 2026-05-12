from __future__ import annotations

from pathlib import Path

from . import __about__
from . import datatypes as _types
from . import reader as _reader
from . import validate as _validate
from .parameters import rcParams

__version__ = __about__.__version__
__all__ = ["Study", "__version__", "rcParams"]


class Study:
    """
    Coordinate one FLiPPR analysis from FragPipe output directories.

    A study reads a Limited Proteolysis (LiP) FragPipe result directory,
    optionally reads a matched Trypsin-only (TrP) directory for protein-level
    normalization, stores one or more control-vs-test comparisons, and runs
    those comparisons into `Result` objects.

    The FragPipe output directory must match the selected acquisition method.
    DDA studies expect `combined_ion.tsv`, `combined_protein.tsv`, and
    `experiment_annotation.tsv`. DIA studies expect `ion.tsv`,
    `dia-quant-output/report.pr_matrix.tsv`,
    `dia-quant-output/report.pg_matrix.tsv`, and
    `experiment_annotation.tsv`.

    Sample names should end with replicate numbers, such as `CTRL_1`,
    `CTRL_2`, `DRUG_1`, and `DRUG_2`. When adding a process, pass only the
    condition prefix (`CTRL`, `DRUG`) and the replicate layout.

    Args:
        lip: FragPipe output directory for the LiP experiment.
        trp: Optional FragPipe output directory for the matched TrP experiment.
            Provide this when calling `add_process()` with TrP normalization
            arguments.
        method: FragPipe acquisition and quantification workflow. Use `"dda"`
            for data-dependent acquisition or `"dia"` for data-independent
            acquisition.

    Attributes:
        lip: Validated LiP output path.
        trp: Validated TrP output path, if provided.
        method: Validated acquisition method.
        processes: Registered comparisons keyed by process ID.
        results: Completed results keyed by process ID after `run()`.

    See Also:
        `samples`: Inspect condition names parsed from FragPipe annotations.
        `add_process`: Register a control-vs-test comparison.
        `run`: Execute all registered comparisons.

    Examples:
        Run a LiP-only DDA analysis:

        >>> import flippr
        >>> study = flippr.Study("fragpipe/LiP_DDA", method="dda")
        >>> study.samples
        {'LiP': {'CTRL', 'DRUG'}}
        >>> study.add_process("drug", lip_ctrl="CTRL", lip_test="DRUG", n_rep=3)
        >>> results = study.run()
        >>> results["drug"].peptide

        Run a DIA analysis with matched TrP normalization:

        >>> study = flippr.Study(
        ...     lip="fragpipe/LiP_DIA",
        ...     trp="fragpipe/TrP_DIA",
        ...     method="dia",
        ... )
        >>> study.add_process(
        ...     "drug_norm",
        ...     lip_ctrl="CTRL",
        ...     lip_test="DRUG",
        ...     n_rep=3,
        ...     trp_ctrl="CTRL",
        ...     trp_test="DRUG",
        ...     trp_n_rep=3,
        ... )
        >>> results = study.run()
        >>> results["drug_norm"].protein_summary
    """

    def __init__(
        self,
        lip: str | Path,
        trp: str | Path | None = None,
        method: str = "dda",
    ) -> None:
        _lip, _trp, _method = _validate._validate_study(lip, trp, method)

        self.lip = _lip
        self.trp = _trp
        self.method = _method
        self.processes: dict[str, _types.Process] = {}
        self.results: dict[str, _types.Result] = {}

    @property
    def samples(self) -> dict[str, set[str]]:
        """
        Condition names found in the FragPipe experiment annotations.

        FLiPPR reads `experiment_annotation.tsv`, removes the final replicate
        suffix from each sample name, and returns the remaining condition names.
        For sample names such as `CTRL_1`, `CTRL_2`, and `DRUG_1`, this returns
        `{"CTRL", "DRUG"}`.

        Returns:
            A dictionary with a `LiP` key and, when a TrP path was provided, a
            `TrP` key. Each value is the set of detected condition names.

        Examples:
            >>> study = flippr.Study("fragpipe/LiP_DDA")
            >>> study.samples
            {'LiP': {'CTRL', 'DRUG'}}
        """

        lip_samples = self._get_samples(self.lip)

        if self.trp is not None:
            trp_samples = self._get_samples(self.trp)

            return {"LiP": lip_samples, "TrP": trp_samples}

        return {"LiP": lip_samples}

    def _get_samples(self, path: Path) -> set[str]:
        """
        Parse condition names from one FragPipe experiment annotation file.
        """

        annot = _reader._read_experiment_annotation(path)

        return {
            "_".join(sample.split("_")[:-1])
            for info in annot.values()
            if (sample := info.get("Sample Name"))
        }

    def add_process(
        self,
        pid: str,
        lip_ctrl: str,
        lip_test: str,
        n_rep: _types.Replicate,
        trp_ctrl: str | None = None,
        trp_test: str | None = None,
        trp_n_rep: _types.Replicate | None = None,
    ) -> None:
        """
        Register one control-vs-test comparison.

        A process defines which LiP sample groups should be compared and, when
        a TrP directory was provided to `Study`, which TrP sample groups should
        be used for protein-level normalization. Registered processes are stored
        in `Study.processes` and are executed by `Study.run()`.

        Replicate arguments may be:

        - `3` for equal control and test replicate counts.
        - `(3, 4)` for different control and test replicate counts.
        - `((1, 3, 4), (1, 2, 5))` for explicit replicate suffixes.

        Args:
            pid: Unique process ID used as the key in `Study.processes` and
                `Study.results`.
            lip_ctrl: LiP control condition prefix, such as `"CTRL"`.
            lip_test: LiP test condition prefix, such as `"DRUG"`.
            n_rep: LiP replicate layout.
            trp_ctrl: TrP control condition prefix. Required only for TrP
                normalization.
            trp_test: TrP test condition prefix. Required only for TrP
                normalization.
            trp_n_rep: TrP replicate layout. Required only for TrP
                normalization.

        Raises:
            ValueError: If only part of the TrP normalization arguments are
                provided.
            TypeError: If the replicate layout is not supported.

        Examples:
            Add a LiP-only comparison:

            >>> study.add_process("drug", "CTRL", "DRUG", 3)

            Add a normalized LiP comparison:

            >>> study.add_process(
            ...     "drug_norm",
            ...     lip_ctrl="CTRL",
            ...     lip_test="DRUG",
            ...     n_rep=3,
            ...     trp_ctrl="CTRL",
            ...     trp_test="DRUG",
            ...     trp_n_rep=3,
            ... )

            Add two comparisons before one run:

            >>> study.add_process("low", "CTRL", "DRUG_LOW", 3)
            >>> study.add_process("high", "CTRL", "DRUG_HIGH", 3)
        """

        self.processes.update(
            {
                pid: _types.Process(
                    rcParams,
                    self.lip,
                    self.trp,
                    self.method,
                    pid,
                    lip_ctrl,
                    lip_test,
                    n_rep,
                    trp_ctrl,
                    trp_test,
                    trp_n_rep,
                )
            }
        )

    def run(self) -> dict[str, _types.Result]:
        """
        Execute all registered processes.

        `run()` evaluates each process in insertion order, stores the output in
        `Study.results`, and returns the same dictionary. Set any global options
        in `flippr.rcParams` before calling this method.

        Returns:
            A dictionary mapping process IDs to `Result` objects. Each result
            exposes Polars dataframes such as `ion`, `peptide`, `cut_site`, and
            `protein_summary`.

        Examples:
            >>> flippr.rcParams["ion.missing_intensity_thresh"] = 1
            >>> results = study.run()
            >>> results["drug"].protein_summary
        """

        self.results = {pid: proc.run() for pid, proc in self.processes.items()}

        return self.results
