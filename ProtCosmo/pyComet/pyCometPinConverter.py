#!/usr/bin/env python3
"""Convert Comet PIN files among text/gzip/parquet formats.

This tool supports:
- .pin
- .pin.gz
- .pin.parquet
- .pin.parquet.gz

Behavior summary:
- PIN-text input (.pin/.pin.gz):
  - normalize Peptide to UniMod style by removing flanks and converting [mass] to [U:id]
  - normalize Proteins to one comma-joined field with --proteins-keep
- Parquet input:
  - keep Peptide text as-is
  - normalize Proteins with --proteins-keep
"""

from __future__ import annotations

import argparse
import gzip
import io
import math
import os
import re
import sys
from typing import Callable, Dict, List, Optional, Sequence, Tuple

sys.path.append(os.path.dirname(__file__))

from utils.pin_converter_utils import (
    WarningSink,
    load_unimod_entries_flexible,
    normalize_pin_peptide_to_unimod,
)


PIN_COLUMNS_NO_PROTEINS: List[str] = [
    "SpecId",
    "Label",
    "ScanNr",
    "ExpMass",
    "CalcMass",
    "lnrSp",
    "deltLCn",
    "deltCn",
    "lnExpect",
    "Xcorr",
    "Sp",
    "IonFrac",
    "Mass",
    "PepLen",
    "Charge1",
    "Charge2",
    "Charge3",
    "Charge4",
    "Charge5",
    "Charge6",
    "enzN",
    "enzC",
    "enzInt",
    "lnNumSP",
    "dM",
    "absdM",
    "Peptide",
]
PIN_COLUMNS: List[str] = PIN_COLUMNS_NO_PROTEINS + ["Proteins"]

PARQUET_SCHEMA: Dict[str, str] = {
    "SpecId": "string",
    "Label": "int8",
    "ScanNr": "int32",
    "ExpMass": "float32",
    "CalcMass": "float32",
    "lnrSp": "float32",
    "deltLCn": "int8",
    "deltCn": "float32",
    "lnExpect": "float32",
    "Xcorr": "float32",
    "Sp": "uint16",
    "IonFrac": "float32",
    "Mass": "float32",
    "PepLen": "int8",
    "Charge1": "int8",
    "Charge2": "int8",
    "Charge3": "int8",
    "Charge4": "int8",
    "Charge5": "int8",
    "Charge6": "int8",
    "enzN": "int8",
    "enzC": "int8",
    "enzInt": "uint8",
    "lnNumSP": "float32",
    "dM": "float32",
    "absdM": "float32",
    "Peptide": "string",
    "Proteins": "string",
}

PARQUET_WRITE_KWARGS: Dict[str, object] = {
    "engine": "pyarrow",
    "compression": "zstd",
    "compression_level": 3,
    "use_dictionary": False,
    "index": False,
}

SUPPORTED_OUTPUT_FORMATS = (
    ".pin",
    ".pin.gz",
    ".pin.parquet",
    ".pin.parquet.gz",
    ".tsv",
    ".tsv.gz",
    ".parquet",
    ".parquet.gz",
)
_INT_TOKEN_RE = re.compile(r"^[+-]?\d+$")
_FLOAT_TOKEN_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_UNIMOD_PATH = os.path.join(_SCRIPT_DIR, "resource", "common_unimod.csv")


class ConversionError(Exception):
    """Raised when conversion cannot proceed."""


def detect_output_format(path: str) -> str:
    lower = path.lower()
    if lower.endswith(".pin.parquet.gz") or lower.endswith(".parquet.gz"):
        return "parquet_gz"
    if lower.endswith(".pin.parquet") or lower.endswith(".parquet"):
        return "parquet"
    if lower.endswith(".pin.gz") or lower.endswith(".tsv.gz"):
        return "pin_gz"
    if lower.endswith(".pin") or lower.endswith(".tsv"):
        return "pin"
    raise ConversionError(
        f"unsupported file extension for '{path}'. supported: {', '.join(SUPPORTED_OUTPUT_FORMATS)}"
    )


def detect_input_format(path: str, input_percolator: bool) -> str:
    if input_percolator:
        lower = path.lower()
        if lower.endswith(".parquet.gz"):
            return "parquet_gz"
        if lower.endswith(".parquet"):
            return "parquet"
        if lower.endswith(".gz"):
            return "percolator_tsv_gz"
        return "percolator_tsv"
    return detect_output_format(path)


def default_output_path(input_path: str, input_format: str) -> str:
    if input_format in ("parquet", "parquet_gz"):
        raise ConversionError(
            "--output is required when input is .parquet or .parquet.gz"
        )
    if input_format in ("pin_gz", "percolator_tsv_gz"):
        # x.<text>.gz -> x.<text>.parquet
        return input_path[: -len(".gz")] + ".parquet"
    # x.<text> -> x.<text>.parquet
    return input_path + ".parquet"


def _open_text_reader(path: str, is_gzip: bool):
    if is_gzip:
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def _open_text_writer(path: str, is_gzip: bool):
    if is_gzip:
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return open(path, "w", encoding="utf-8", newline="")


def _normalize_proteins(parts: Sequence[object], proteins_keep: int) -> str:
    values: List[str] = []
    for part in parts:
        if part is None:
            continue
        text = str(part).strip()
        if not text:
            continue
        split_items = re.split(r"[,;]", text)
        values.extend(item.strip() for item in split_items if item.strip())

    if proteins_keep > 0:
        values = values[:proteins_keep]
    return ",".join(values)


def _find_column_case_insensitive(
    columns: Sequence[str],
    candidate_names: Sequence[str],
) -> Optional[str]:
    lookup: Dict[str, str] = {}
    for col in columns:
        lookup[str(col).lower()] = str(col)
    for name in candidate_names:
        matched = lookup.get(name.lower())
        if matched is not None:
            return matched
    return None


def read_pin_text_rows(
    path: str,
    is_gzip: bool,
    proteins_keep: int,
    peptide_mapper: Callable[[str], str],
) -> Tuple[List[str], List[Dict[str, str]]]:
    rows: List[Dict[str, str]] = []
    with _open_text_reader(path, is_gzip) as handle:
        header_line = handle.readline()
        if not header_line:
            raise ConversionError(f"input file is empty: {path}")
        header = header_line.rstrip("\r\n").split("\t")
        if "Peptide" not in header:
            raise ConversionError(f"invalid PIN header in {path}: missing required 'Peptide' column")
        if len(header) >= len(PIN_COLUMNS_NO_PROTEINS) and (
            header[: len(PIN_COLUMNS_NO_PROTEINS)] != PIN_COLUMNS_NO_PROTEINS
        ):
            print(
                "warning: PIN header does not match expected first 27 columns; parsing by position.",
                file=sys.stderr,
            )
        proteins_idx = header.index("Proteins") if "Proteins" in header else -1
        output_columns = list(header)
        if proteins_idx < 0:
            output_columns.append("Proteins")
        warned_proteins_not_last = False

        for line in handle:
            stripped = line.rstrip("\r\n")
            if not stripped:
                continue
            parts = stripped.split("\t")
            row: Dict[str, str] = {}
            for idx, col in enumerate(header):
                row[col] = parts[idx] if idx < len(parts) else ""
            trailing = parts[len(header) :] if len(parts) > len(header) else []

            proteins_raw: List[str] = []
            if proteins_idx >= 0 and proteins_idx < len(parts):
                proteins_raw.append(parts[proteins_idx])
            if trailing:
                proteins_raw.extend(trailing)
            if proteins_idx >= 0 and proteins_idx != len(header) - 1 and trailing:
                if not warned_proteins_not_last:
                    print(
                        "warning: 'Proteins' is not the last header column but row has extra trailing fields; "
                        "trailing fields are appended into Proteins.",
                        file=sys.stderr,
                    )
                    warned_proteins_not_last = True
            row["Proteins"] = _normalize_proteins(proteins_raw, proteins_keep)
            row["Peptide"] = peptide_mapper(row.get("Peptide", ""))
            rows.append(row)
    return output_columns, rows


def read_percolator_text_rows(
    path: str,
    is_gzip: bool,
    proteins_keep: int,
    peptide_mapper: Callable[[str], str],
) -> Tuple[List[str], List[Dict[str, str]]]:
    rows: List[Dict[str, str]] = []
    with _open_text_reader(path, is_gzip) as handle:
        header_line = handle.readline()
        if not header_line:
            raise ConversionError(f"input file is empty: {path}")
        header = header_line.rstrip("\r\n").split("\t")
        peptide_col = _find_column_case_insensitive(header, ("peptide", "Peptide"))
        if peptide_col is None:
            raise ConversionError(
                f"invalid percolator header in {path}: missing required 'peptide' column"
            )
        protein_col = _find_column_case_insensitive(
            header,
            ("proteinIds", "proteinids", "Proteins", "proteins"),
        )
        if protein_col is None:
            raise ConversionError(
                f"invalid percolator header in {path}: missing required 'proteinIds' column"
            )
        output_columns = list(header)

        warned_protein_not_last = False
        protein_idx = header.index(protein_col)

        for line in handle:
            stripped = line.rstrip("\r\n")
            if not stripped:
                continue
            parts = stripped.split("\t")
            row: Dict[str, str] = {}
            for idx, col in enumerate(header):
                row[col] = parts[idx] if idx < len(parts) else ""
            trailing = parts[len(header) :] if len(parts) > len(header) else []

            proteins_raw: List[str] = []
            proteins_raw.append(row.get(protein_col, ""))
            if trailing:
                proteins_raw.extend(trailing)
            if protein_idx != len(header) - 1 and trailing:
                if not warned_protein_not_last:
                    print(
                        "warning: 'proteinIds' is not the last header column but row has extra trailing fields; "
                        "trailing fields are appended into proteinIds.",
                        file=sys.stderr,
                    )
                    warned_protein_not_last = True
            normalized_proteins = _normalize_proteins(proteins_raw, proteins_keep)
            row[protein_col] = normalized_proteins

            row[peptide_col] = peptide_mapper(row.get(peptide_col, ""))
            rows.append(row)
    return output_columns, rows


def write_pin_text_rows(
    path: str,
    is_gzip: bool,
    rows: Sequence[Dict[str, object]],
    columns: Sequence[str],
) -> None:
    with _open_text_writer(path, is_gzip) as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            fields: List[str] = []
            for col in columns:
                value = row.get(col, "")
                fields.append(_format_pin_scalar(value))
            handle.write("\t".join(fields) + "\n")


def _is_missing_value(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    # pandas.NA without hard dependency.
    if type(value).__name__ == "NAType":
        return True
    return False


def _format_pin_scalar(value: object) -> str:
    if _is_missing_value(value):
        return ""
    if isinstance(value, bool):
        return "1" if value else "0"
    return str(value)


def _ensure_parquet_dependencies():
    try:
        import pandas as pd  # type: ignore
        import pyarrow  # noqa: F401
    except ImportError as exc:
        raise ConversionError(
            "parquet conversion requires pandas and pyarrow. "
            "please install them in this environment."
        ) from exc
    return pd


def _is_blank_like(value: object) -> bool:
    if _is_missing_value(value):
        return True
    if isinstance(value, str) and not value.strip():
        return True
    return False


def _infer_unknown_dtype(series) -> str:
    dtype_name = str(series.dtype)
    if dtype_name.startswith("int") or dtype_name.startswith("uint") or dtype_name.startswith("Int"):
        return "int32"
    if dtype_name.startswith("float"):
        return "float32"

    values: List[str] = []
    for value in series.tolist():
        if _is_blank_like(value):
            continue
        values.append(str(value).strip())
    if not values:
        return "string"
    if all(_INT_TOKEN_RE.match(text) for text in values):
        return "int32"
    if all(_FLOAT_TOKEN_RE.match(text) for text in values):
        return "float32"
    return "string"


def _cast_series_to_dtype(series, dtype: str):
    pd = _ensure_parquet_dependencies()
    if dtype == "string":
        return series.astype("string")

    as_text = series.map(lambda value: "" if _is_blank_like(value) else str(value).strip())
    numeric = pd.to_numeric(series, errors="coerce")
    invalid_mask = (as_text != "") & numeric.isna()
    if bool(invalid_mask.any()):
        raise ConversionError(
            "cannot cast column to numeric dtype due to non-numeric values"
        )

    if dtype == "float32":
        return numeric.astype("float32")
    if dtype in ("int8", "uint8", "uint16", "int32"):
        if numeric.isna().any():
            if dtype == "int32":
                return numeric.astype("Int32")
            raise ConversionError(f"column has blank values that cannot be cast to {dtype}")
        return numeric.astype(dtype)
    raise ConversionError(f"unsupported dtype cast request: {dtype}")


def _cast_dataframe_to_output_schema(df, column_order: Sequence[str]):
    converted = df.copy()
    for col in column_order:
        if col not in converted.columns:
            converted[col] = ""
    converted = converted[list(column_order)]

    for col in column_order:
        target_dtype = PARQUET_SCHEMA.get(col)
        if target_dtype is None:
            target_dtype = _infer_unknown_dtype(converted[col])
        converted[col] = _cast_series_to_dtype(converted[col], target_dtype)
    return converted


def _rows_to_dataframe(rows: Sequence[Dict[str, object]], columns: Sequence[str]):
    pd = _ensure_parquet_dependencies()
    df = pd.DataFrame(list(rows), columns=list(columns))
    return _cast_dataframe_to_output_schema(df, columns)


def _read_parquet_dataframe(path: str, file_format: str):
    pd = _ensure_parquet_dependencies()
    if file_format == "parquet":
        df = pd.read_parquet(path, engine="pyarrow")
    else:
        with gzip.open(path, "rb") as handle:
            payload = handle.read()
        df = pd.read_parquet(io.BytesIO(payload), engine="pyarrow")
    df.columns = [str(col) for col in df.columns]
    return df


def _write_parquet_dataframe(df, path: str, file_format: str, column_order: Sequence[str]) -> None:
    casted = _cast_dataframe_to_output_schema(df, column_order)
    if file_format == "parquet":
        casted.to_parquet(path, **PARQUET_WRITE_KWARGS)
        return

    buffer = io.BytesIO()
    casted.to_parquet(buffer, **PARQUET_WRITE_KWARGS)
    with gzip.open(path, "wb") as handle:
        handle.write(buffer.getvalue())


def _normalize_dataframe_proteins(df, proteins_keep: int):
    normalized = df.copy()
    if "Proteins" not in normalized.columns:
        normalized["Proteins"] = ""
    normalized["Proteins"] = normalized["Proteins"].map(
        lambda value: _normalize_proteins([value], proteins_keep)
    )
    return normalized


def _dataframe_to_pin_rows(df, columns: Sequence[str]) -> List[Dict[str, object]]:
    prepared = df.copy()
    for col in columns:
        if col not in prepared.columns:
            prepared[col] = ""
    prepared = prepared[list(columns)]
    records = prepared.to_dict(orient="records")
    rows: List[Dict[str, object]] = []
    for record in records:
        row: Dict[str, object] = {}
        for col in columns:
            row[col] = record.get(col, "")
        rows.append(row)
    return rows


def _load_unimod_entries_required(unimod_path: str):
    if not os.path.isfile(unimod_path):
        raise ConversionError(f"unimod file does not exist: {unimod_path}")
    entries = load_unimod_entries_flexible(unimod_path)
    if not entries:
        raise ConversionError(
            f"cannot load unimod entries from: {unimod_path}. expected first 3 columns: mass,residues,unimod_id"
        )
    return entries


def _normalize_percolator_dataframe(
    df,
    proteins_keep: int,
    peptide_mapper: Callable[[str], str],
):
    normalized = df.copy()
    columns = [str(col) for col in normalized.columns]
    peptide_col = _find_column_case_insensitive(columns, ("peptide", "Peptide"))
    if peptide_col is None:
        raise ConversionError(
            "invalid percolator parquet input: missing required 'peptide' column"
        )
    protein_col = _find_column_case_insensitive(
        columns,
        ("proteinIds", "proteinids", "Proteins", "proteins"),
    )
    if protein_col is None:
        raise ConversionError(
            "invalid percolator parquet input: missing required 'proteinIds' column"
        )

    normalized[peptide_col] = normalized[peptide_col].map(
        lambda value: peptide_mapper("" if _is_blank_like(value) else str(value))
    )

    normalized[protein_col] = normalized[protein_col].map(
        lambda value: _normalize_proteins([value], proteins_keep)
    )
    return normalized


def convert_peptide_string(
    peptide_text: str,
    unimod_path: str = _DEFAULT_UNIMOD_PATH,
    unimod_ppm: float = 10.0,
    warning_sink: Optional[WarningSink] = None,
) -> str:
    """Convert one PIN-style peptide string to UniMod style.

    Example:
    - input:  -.n[42.0106]M[15.9949]LQFLLEVNK.S
    - output: [U:1]-M[U:35]LQFLLEVNK
    """
    if unimod_ppm < 0.0:
        raise ConversionError("--unimod-ppm must be >= 0")
    entries = _load_unimod_entries_required(unimod_path)
    sink = warning_sink if warning_sink is not None else WarningSink()
    return normalize_pin_peptide_to_unimod(peptide_text, entries, unimod_ppm, sink)


def convert_file(
    input_path: str,
    output_path: str,
    proteins_keep: int,
    unimod_path: str,
    unimod_ppm: float,
    input_percolator: bool = False,
) -> None:
    input_format = detect_input_format(input_path, input_percolator)
    output_format = detect_output_format(output_path)

    if proteins_keep < 0:
        raise ConversionError("--proteins-keep must be >= 0")
    if unimod_ppm < 0.0:
        raise ConversionError("--unimod-ppm must be >= 0")

    if input_percolator:
        unimod_entries = _load_unimod_entries_required(unimod_path)
        warning_sink = WarningSink()
        peptide_mapper = lambda peptide: normalize_pin_peptide_to_unimod(
            peptide,
            unimod_entries,
            unimod_ppm,
            warning_sink,
        )

        if input_format in ("percolator_tsv", "percolator_tsv_gz"):
            columns, rows = read_percolator_text_rows(
                path=input_path,
                is_gzip=(input_format == "percolator_tsv_gz"),
                proteins_keep=proteins_keep,
                peptide_mapper=peptide_mapper,
            )
            if output_format in ("pin", "pin_gz"):
                write_pin_text_rows(
                    output_path,
                    output_format == "pin_gz",
                    rows,
                    columns,
                )
                return
            df = _rows_to_dataframe(rows, columns)
            _write_parquet_dataframe(df, output_path, output_format, columns)
            return

        # percolator parquet input
        df_in = _read_parquet_dataframe(input_path, input_format)
        df_in = _normalize_percolator_dataframe(
            df_in,
            proteins_keep=proteins_keep,
            peptide_mapper=peptide_mapper,
        )
        columns = [str(col) for col in df_in.columns]
        if output_format in ("parquet", "parquet_gz"):
            _write_parquet_dataframe(df_in, output_path, output_format, columns)
            return
        rows_out = _dataframe_to_pin_rows(df_in, columns)
        write_pin_text_rows(output_path, output_format == "pin_gz", rows_out, columns)
        return

    if input_format in ("pin", "pin_gz"):
        unimod_entries = _load_unimod_entries_required(unimod_path)

        warning_sink = WarningSink()
        columns, rows = read_pin_text_rows(
            path=input_path,
            is_gzip=(input_format == "pin_gz"),
            proteins_keep=proteins_keep,
            peptide_mapper=lambda peptide: normalize_pin_peptide_to_unimod(
                peptide,
                unimod_entries,
                unimod_ppm,
                warning_sink,
            ),
        )

        if output_format in ("pin", "pin_gz"):
            write_pin_text_rows(
                output_path,
                output_format == "pin_gz",
                rows,
                columns,
            )
            return

        df = _rows_to_dataframe(rows, columns)
        _write_parquet_dataframe(df, output_path, output_format, columns)
        return

    # parquet input
    df_in = _read_parquet_dataframe(input_path, input_format)
    df_in = _normalize_dataframe_proteins(df_in, proteins_keep)
    columns = [str(col) for col in df_in.columns]

    if output_format in ("parquet", "parquet_gz"):
        _write_parquet_dataframe(df_in, output_path, output_format, columns)
        return

    # parquet -> pin text, keep peptide strings as-is
    rows_out = _dataframe_to_pin_rows(df_in, columns)
    write_pin_text_rows(output_path, output_format == "pin_gz", rows_out, columns)


class CometHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def _normalize_argv_for_input_str(argv: Sequence[str]) -> List[str]:
    """Allow --input-str values that start with '-.' without requiring '=' syntax."""
    out = list(argv)
    idx = 0
    while idx < len(out):
        token = out[idx]
        if token != "--input-str":
            idx += 1
            continue
        if idx + 1 >= len(out):
            return out
        candidate = out[idx + 1]
        if candidate.startswith("-."):
            out[idx] = f"--input-str={candidate}"
            del out[idx + 1]
        idx += 1
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Comet PIN files among .pin, .pin.gz, .pin.parquet, .pin.parquet.gz.\n"
            "For PIN-text input, Peptide is normalized to UniMod tags ([U:id]) and flanks are removed.\n"
            "For parquet input, Peptide text is preserved.\n"
            "With --input-percolator, non-.parquet input is treated as percolator TSV."
        ),
        epilog=(
            "Custom --unimod file format (CSV or TSV):\n"
            "  first 3 columns must be: mass_delta, residues, unimod_id\n"
            "  extra columns are optional names/notes\n"
            "  example row:\n"
            "    +15.994915,DKNPFYRMCWHGUEILQSTV,UNIMOD:35,Oxidation\n\n"
            "Examples:\n"
            "  python pyCometPinConverter.py --input x.pin\n"
            "  python pyCometPinConverter.py --input x.pin.gz --output x.pin.parquet.gz\n"
            "  python pyCometPinConverter.py --input x.pin.parquet --output x.pin\n"
            "  python pyCometPinConverter.py --input xx.decoy.psms.tsv --input-percolator\n"
            "  python pyCometPinConverter.py --input-str '-.n[42.0106]M[15.9949]LQFLLEVNK.S'\n"
            "  python pyCometPinConverter.py --input x.pin --unimod custom_unimod.tsv --unimod-ppm 10\n"
        ),
        formatter_class=CometHelpFormatter,
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input", help="Input file path.")
    input_group.add_argument(
        "--input-str",
        help="Single PIN-style peptide string to convert and print.",
    )
    parser.add_argument(
        "--input-percolator",
        action="store_true",
        help=(
            "Treat --input as percolator output. "
            "If input path does not end with .parquet/.parquet.gz, it is parsed as TSV. "
            "Column names are preserved; only peptide/proteinIds values are normalized."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Output file path. If omitted and input is text (.pin/.pin.gz or --input-percolator TSV), "
            "output defaults to <input>.parquet. For parquet input, --output is required."
        ),
    )
    parser.add_argument(
        "--proteins-keep",
        type=int,
        default=1,
        help="Keep first N proteins in Proteins column and join with comma. 0 keeps all.",
    )
    parser.add_argument(
        "--unimod",
        default=_DEFAULT_UNIMOD_PATH,
        help="UniMod CSV/TSV path. Expected first 3 columns: mass_delta, residues, unimod_id.",
    )
    parser.add_argument(
        "--unimod-ppm",
        type=float,
        default=10.0,
        help="PPM tolerance for mapping [mass] to [U:id].",
    )
    argv = _normalize_argv_for_input_str(sys.argv[1:])
    args = parser.parse_args(argv)

    if args.proteins_keep < 0:
        parser.error("--proteins-keep must be >= 0")
    if args.unimod_ppm < 0.0:
        parser.error("--unimod-ppm must be >= 0")

    if args.input_str is not None:
        if args.input_percolator:
            parser.error("--input-percolator cannot be used with --input-str")
        if args.output:
            parser.error("--output cannot be used with --input-str")
        try:
            converted = convert_peptide_string(
                peptide_text=args.input_str,
                unimod_path=args.unimod,
                unimod_ppm=args.unimod_ppm,
            )
        except ConversionError as exc:
            print(f"error: {exc}", file=sys.stderr)
            raise SystemExit(2) from exc
        print(converted)
        return

    input_format = detect_input_format(args.input, args.input_percolator)
    output_path = args.output if args.output else default_output_path(args.input, input_format)

    try:
        convert_file(
            input_path=args.input,
            output_path=output_path,
            proteins_keep=args.proteins_keep,
            unimod_path=args.unimod,
            unimod_ppm=args.unimod_ppm,
            input_percolator=args.input_percolator,
        )
    except ConversionError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
