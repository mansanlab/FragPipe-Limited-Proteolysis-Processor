from __future__ import annotations

from pathlib import Path
from typing import Literal, cast
from warnings import warn

from .parameters import _DDA_FP_FILES, _DIA_FP_FILES

type Method = Literal["dda", "dia"]
type Replicate = int | tuple[int, int] | tuple[tuple[int, ...], tuple[int, ...]]
type ReplicateKind = Literal["int", "tuple", "tuple_tuple"]


def _validate_study(
    lip: str | Path,
    trp: str | Path | None,
    method: str,
) -> tuple[Path, Path | None, Method]:
    method = __validate_method(method)

    lip = __validate_fragpipe_path(lip, "lip", method)

    if trp is not None:
        trp = __validate_fragpipe_path(trp, "trp", method)

    return lip, trp, method


def __validate_method(method: str) -> Method:
    if not isinstance(method, str):
        raise TypeError(
            f'`method` was provided: "{method}" with type `{type(method)}`. '
            'Set `method` to "dda" or "dia".'
        )

    if method not in {"dda", "dia"}:
        raise ValueError(
            f'`method` was provided: "{method}". "{method}" is not recognized. '
            'Set `method` to "dda" or "dia".'
        )

    return cast(Method, method)


def __validate_fragpipe_path(path: str | Path, liptrp: str, method: Method) -> Path:
    if not isinstance(path, (str, Path)):
        raise TypeError(
            f'`{liptrp}` was provided: "{path}" with type `{type(path)}`. '
            f"Set `{liptrp}` to a FragPipe output path."
        )

    path = Path(path).expanduser()

    if not path.is_dir():
        raise ValueError(
            f'`{liptrp}` was provided: "{path}". "{path}" is not a directory path. '
            f"Set `{liptrp}` to a FragPipe output directory path."
        )

    __validate_fragpipe_files(path, method)

    return path


def __validate_fragpipe_files(path: Path, method: Method) -> None:
    required = _DDA_FP_FILES if method == "dda" else _DIA_FP_FILES

    if all((path / file).exists() for file in required):
        return

    missing = [file for file in required if not (path / file).exists()]
    files = "\n".join(f"- `{file}`" for file in missing)
    raise FileNotFoundError(
        f'Files not found in "{path}". Missing required FragPipe output files:\n{files}'
    )


def _validate_replicate(replicate: Replicate) -> ReplicateKind:
    def __int_validation(i: int) -> None:
        if isinstance(i, bool):
            raise TypeError("Replicate values must be integers, not booleans.")

        if i <= 1:
            raise ValueError(
                f"Number of replicates was set to `{i}`. "
                "Replicate values must be positive and greater than or equal to `2`."
            )

        if i == 2:
            warn(
                "Using two replicates will result in an under-powered study.",
                UserWarning,
                stacklevel=2,
            )
            warn(
                "Set flippr.rcParams `ion.missing_intensity_thresh` to `0` for the best results.",
                UserWarning,
                stacklevel=2,
            )

    if isinstance(replicate, int):
        __int_validation(replicate)
        return "int"

    if isinstance(replicate, tuple):
        replicate = cast(tuple, replicate)

        if len(replicate) != 2:
            raise ValueError(
                f"Replicate input contains `{len(replicate)}` `{type(replicate)}`. "
                "Only `int`, `tuple[int, int]`, or "
                "`tuple[tuple[int, ...], tuple[int, ...]]` are allowed."
            )

        if all(isinstance(i, int) for i in replicate):
            replicate = cast(tuple[int, int], replicate)

            for value in replicate:
                __int_validation(value)
            return "tuple"

        if all(isinstance(i, tuple) for i in replicate):
            replicate = cast(tuple[tuple[int, ...], tuple[int, ...]], replicate)

            if all(
                all(isinstance(i, int) and not isinstance(i, bool) for i in tup)
                for tup in replicate
            ):
                for tup in replicate:
                    __int_validation(len(tup))
                return "tuple_tuple"

    raise TypeError(
        f"Replicate input is of type `{type(replicate)}`. "
        "Only `int`, `tuple[int, int]`, or "
        "`tuple[tuple[int, ...], tuple[int, ...]]` are allowed."
    )
