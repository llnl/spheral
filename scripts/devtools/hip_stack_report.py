#!/usr/bin/env python3
"""
Summarize HIP kernel stack metadata from ROCm saved-temp assembly.

Build with:

  -DSPHERAL_HIP_REPORT_STACK_USAGE=ON

Then run this script on the CMake build directory.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, List, Optional, Sequence, TextIO


METADATA_SUFFIXES = {".s", ".yaml", ".yml"}
SUMMARY_HEADERS = [
    "component",
    "kernels",
    "variable-stack kernels",
    "largest fixed stack",
    "largest fixed stack in variable kernels",
]

KERNELS_RE = re.compile(r"^(?P<indent>\s*)amdhsa\.kernels:\s*$")
NAME_RE = re.compile(r"^(?P<indent>\s*)(?P<list_item>-\s+)?\.name:\s*(?P<value>.+?)\s*$")
FIELD_RE = re.compile(r"^(?P<indent>\s*)(?P<key>\.[A-Za-z0-9_]+):\s*(?P<value>.+?)\s*$")
SAVED_TEMP_SUFFIX_RE = re.compile(r"-hip-amdgcn-amd-amdhsa(?:-[A-Za-z0-9_.+]+)?$")


@dataclass
class KernelRecord:
    source: str
    line: int
    kernel: str
    static_b: Optional[int] = None
    dynamic: Optional[bool] = None
    sgprs: Optional[int] = None
    vgprs: Optional[int] = None

    def has_stack_metadata(self) -> bool:
        return self.static_b is not None or self.dynamic is not None

    def sort_key(self) -> Any:
        return (
            self.static_b if self.static_b is not None else -1,
            self.dynamic is True,
            self.sgprs if self.sgprs is not None else -1,
            self.vgprs if self.vgprs is not None else -1,
            self.source,
            self.line,
            self.kernel,
        )


@dataclass
class InputSummary:
    label: str
    records: List[KernelRecord]


def clean_value(value: str) -> str:
    value = value.strip().split("#", 1)[0].strip().rstrip(",")
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]
    return value


def clean_kernel_name(value: str) -> str:
    value = clean_value(value)
    if value.startswith("*Z"):
        return "_" + value[1:]
    return value


def parse_int(value: str) -> Optional[int]:
    value = clean_value(value)
    try:
        return int(value, 0)
    except ValueError:
        match = re.search(r"-?\d+", value)
        return int(match.group(0), 10) if match else None


def parse_bool(value: str) -> Optional[bool]:
    value = clean_value(value).lower()
    if value in ("true", "yes", "1"):
        return True
    if value in ("false", "no", "0"):
        return False
    return None


def device_assembly(path: Path) -> bool:
    return path.suffix.lower() == ".s" and "hip-amdgcn-amd-amdhsa" in path.name


def input_files(paths: Iterable[Path]) -> Iterator[Path]:
    seen: set[str] = set()

    def yield_once(path: Path) -> Iterator[Path]:
        key = str(path.resolve())
        if key not in seen:
            seen.add(key)
            yield path

    for path in paths:
        if path.is_dir():
            for child in path.rglob("*"):
                if child.is_file() and (device_assembly(child) or child.suffix.lower() in (".yaml", ".yml")):
                    yield from yield_once(child)
        elif path.is_file():
            if path.suffix.lower() in METADATA_SUFFIXES:
                yield from yield_once(path)
        else:
            print(f"warning: input path does not exist: {path}", file=sys.stderr)


def append_record(records: List[KernelRecord], record: Optional[KernelRecord]) -> None:
    if record is not None and record.has_stack_metadata():
        records.append(record)


def apply_field(record: KernelRecord, key: str, value: str) -> None:
    key = key.lstrip(".")
    if key == "private_segment_fixed_size":
        record.static_b = parse_int(value)
    elif key == "uses_dynamic_stack":
        record.dynamic = parse_bool(value)
    elif key == "sgpr_count":
        record.sgprs = parse_int(value)
    elif key == "vgpr_count":
        record.vgprs = parse_int(value)


def parse_metadata(path: Path) -> List[KernelRecord]:
    records: List[KernelRecord] = []
    current: Optional[KernelRecord] = None
    in_kernels = False
    kernels_indent: Optional[int] = None
    kernel_indent: Optional[int] = None

    try:
        lines = path.open(encoding="utf-8", errors="ignore")
    except OSError as exc:
        print(f"warning: cannot read {path}: {exc}", file=sys.stderr)
        return records

    with lines:
        for lineno, line in enumerate(lines, 1):
            kernels_match = KERNELS_RE.match(line)
            if kernels_match:
                append_record(records, current)
                current = None
                in_kernels = True
                kernels_indent = len(kernels_match.group("indent"))
                kernel_indent = None
                continue

            if not in_kernels:
                continue

            stripped = line.strip()
            line_indent = len(line) - len(line.lstrip())
            if (
                stripped
                and kernels_indent is not None
                and line_indent <= kernels_indent
                and not stripped.startswith(("-", "."))
            ):
                append_record(records, current)
                current = None
                in_kernels = False
                kernels_indent = None
                kernel_indent = None
                continue

            name_match = NAME_RE.match(line)
            if name_match:
                append_record(records, current)
                kernel_indent = len(name_match.group("indent"))
                if name_match.group("list_item"):
                    kernel_indent += 2
                current = KernelRecord(
                    source=str(path),
                    line=lineno,
                    kernel=clean_kernel_name(name_match.group("value")),
                )
                continue

            if current is None or kernel_indent is None:
                continue

            field_match = FIELD_RE.match(line)
            if field_match and len(field_match.group("indent")) == kernel_indent:
                apply_field(current, field_match.group("key"), field_match.group("value"))

    append_record(records, current)
    return records


def demangle(records: Sequence[KernelRecord]) -> None:
    cxxfilt = shutil.which("llvm-cxxfilt") or shutil.which("c++filt")
    if cxxfilt is None or not records:
        return

    kernel_names = list(dict.fromkeys(record.kernel for record in records))

    try:
        result = subprocess.run(
            [cxxfilt],
            check=False,
            input="\n".join(kernel_names),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except OSError:
        return

    demangled = {
        kernel_name: demangled_name
        for kernel_name, demangled_name in zip(kernel_names, result.stdout.splitlines())
        if demangled_name
    }
    for record in records:
        record.kernel = demangled.get(record.kernel, record.kernel)


def csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    return value


def byte_value(value: Optional[int]) -> str:
    return "" if value is None else f"{value} B"


def write_tsv(records: Sequence[KernelRecord], output: TextIO = sys.stdout) -> None:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(["source", "line", "static_B", "dynamic", "sgprs", "vgprs", "kernel"])
    for record in records:
        writer.writerow(
            [
                record.source,
                record.line,
                csv_value(record.static_b),
                csv_value(record.dynamic),
                csv_value(record.sgprs),
                csv_value(record.vgprs),
                record.kernel,
            ]
        )


def path_after_source_anchor(parts: Sequence[str], stop: int) -> Sequence[str]:
    anchors = {"src", "tests"}
    starts = [index for index, part in enumerate(parts[:stop]) if part in anchors]
    if starts:
        return parts[max(starts):stop]
    return parts[max(0, stop - 2):stop]


def saved_temp_source_name(path: Path) -> str:
    return SAVED_TEMP_SUFFIX_RE.sub("", path.stem)


def join_label(parts: Sequence[str]) -> str:
    return "/".join(part for part in parts if part)


def cmake_target_dir(path: Path) -> Optional[Path]:
    parts = path.parts
    for cmake_files_index, part in enumerate(parts[:-1]):
        if part == "CMakeFiles" and parts[cmake_files_index + 1].endswith(".dir"):
            return Path(*parts[: cmake_files_index + 2])
    return None


def cmake_target_label_from_dir(path: Path) -> Optional[str]:
    parts = path.parts
    try:
        cmake_files_index = parts.index("CMakeFiles")
    except ValueError:
        return None

    source_parts = path_after_source_anchor(parts, cmake_files_index)
    target = parts[cmake_files_index + 1]
    if target.endswith(".dir"):
        target = target[:-4]
    return join_label([*source_parts, target]) or None


def cmake_target_label(path: Path) -> Optional[str]:
    target_dir = cmake_target_dir(path)
    if target_dir is None:
        return None
    return cmake_target_label_from_dir(target_dir)


def input_label(path: Path) -> str:
    label = cmake_target_label(path)
    if label is not None:
        if path.suffix.lower() in METADATA_SUFFIXES:
            return join_label([label, saved_temp_source_name(path)])
        return label

    if path.suffix.lower() in METADATA_SUFFIXES:
        return saved_temp_source_name(path)

    return path.name or str(path)


def is_cmake_target_dir(path: Path) -> bool:
    return path == cmake_target_dir(path)


def find_cmake_build_root(path: Path) -> Optional[Path]:
    start = path if path.is_dir() else path.parent
    for candidate in (start, *start.parents):
        if (candidate / "CMakeFiles" / "TargetDirectories.txt").is_file():
            return candidate
    return None


def compiled_cmake_target_dir(path: Path) -> bool:
    return (
        path.name.endswith(".dir")
        and (
            (path / "DependInfo.cmake").is_file()
            or (path / "link.txt").is_file()
            or any(path.glob("*DependInfo.json"))
        )
    )


def cmake_target_labels(path: Path) -> dict[str, str]:
    build_root = find_cmake_build_root(path)
    if build_root is None:
        return {}

    target_directories = build_root / "CMakeFiles" / "TargetDirectories.txt"
    labels: dict[str, str] = {}

    try:
        lines = target_directories.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError as exc:
        print(f"warning: cannot read {target_directories}: {exc}", file=sys.stderr)
        return labels

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        target_dir = Path(line)
        if not target_dir.is_absolute():
            target_dir = build_root / target_dir
        if not compiled_cmake_target_dir(target_dir):
            continue

        label = cmake_target_label_from_dir(target_dir)
        if label is not None:
            labels[str(target_dir.resolve())] = label

    return labels


def target_label_for_metadata(path: Path, target_labels: dict[str, str]) -> Optional[str]:
    target_dir = cmake_target_dir(path)
    if target_dir is None:
        return None

    if target_labels:
        return target_labels.get(str(target_dir.resolve()))

    return cmake_target_label_from_dir(target_dir)


def sort_records(records: List[KernelRecord]) -> List[KernelRecord]:
    records.sort(key=lambda record: record.sort_key(), reverse=True)
    return records


def collect_records(paths: Iterable[Path]) -> List[KernelRecord]:
    records: List[KernelRecord] = []
    for path in input_files(paths):
        records.extend(parse_metadata(path))
    return sort_records(records)


def collect_summary(path: Path) -> List[InputSummary]:
    if path.is_file() or is_cmake_target_dir(path):
        records = collect_records([path])
        return [InputSummary(input_label(path), records)] if records else []

    target_labels = cmake_target_labels(path)
    groups: dict[str, List[KernelRecord]] = {}
    for metadata_file in input_files([path]):
        label = target_label_for_metadata(metadata_file, target_labels)
        if label is None:
            continue
        records = parse_metadata(metadata_file)
        if not records:
            continue
        groups.setdefault(label, []).extend(records)

    return [InputSummary(label, sort_records(records)) for label, records in sorted(groups.items())]


def collect_summaries(paths: Iterable[Path]) -> List[InputSummary]:
    summaries: List[InputSummary] = []
    for path in paths:
        summaries.extend(collect_summary(path))
    return summaries


def summary_row(label: str, records: Sequence[KernelRecord]) -> List[Any]:
    fixed_records = [record for record in records if record.static_b is not None]
    dynamic_records = [record for record in records if record.dynamic is True]
    max_static_b = max((record.static_b for record in fixed_records), default=None)
    max_dynamic_static_b = max(
        (record.static_b for record in dynamic_records if record.static_b is not None),
        default=None,
    )
    return [
        label,
        len(records),
        len(dynamic_records),
        byte_value(max_static_b),
        byte_value(max_dynamic_static_b),
    ]


def format_summary_row(row: Sequence[Any], widths: Sequence[int]) -> str:
    cells = [str(value) for value in row]
    formatted = [cells[0].ljust(widths[0])]
    formatted.extend(cell.rjust(width) for cell, width in zip(cells[1:], widths[1:]))
    return "  ".join(formatted)


def write_summary_table(summaries: Sequence[InputSummary], output: TextIO = sys.stdout) -> None:
    rows = [summary_row(summary.label, summary.records) for summary in summaries]
    has_dynamic_stack = any(record.dynamic is True for summary in summaries for record in summary.records)

    if len(summaries) > 1:
        records = [record for summary in summaries for record in summary.records]
        rows.append(summary_row("OVERALL", records))

    table_rows = [SUMMARY_HEADERS] + rows
    widths = [max(len(str(row[index])) for row in table_rows) for index in range(len(SUMMARY_HEADERS))]
    separator = ["-" * width for width in widths]

    print(format_summary_row(SUMMARY_HEADERS, widths), file=output)
    print(format_summary_row(separator, widths), file=output)
    for row in rows:
        print(format_summary_row(row, widths), file=output)

    print(file=output)
    print("Fixed-stack columns are compiler-reported known per-thread stack use.", file=output)
    if has_dynamic_stack:
        print(
            "Variable-stack kernels may use additional runtime stack that is not bounded by this metadata.",
            file=output,
        )


def parse_args(argv: Optional[Sequence[str]]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize HIP kernel stack metadata from saved temp files.")
    parser.add_argument("paths", nargs="+", type=Path, help="Saved-temp build directories or metadata files")
    parser.add_argument("--details", action="store_true", help="Write detailed per-record TSV")
    parser.add_argument("--top", type=int, default=None, help="Only include the top N records in detailed TSV")
    args = parser.parse_args(argv)
    if args.top is not None and not args.details:
        parser.error("--top requires --details")
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    if args.details:
        records = collect_records(args.paths)
        if not records:
            print("No HIP kernel stack metadata records found.")
            return 1
        selected_records = records[: args.top] if args.top is not None else records
        demangle(selected_records)
        write_tsv(selected_records)
    else:
        summaries = collect_summaries(args.paths)
        if not summaries:
            print("No HIP kernel stack metadata records found.")
            return 1
        write_summary_table(summaries)

    return 0


if __name__ == "__main__":
    sys.exit(main())
