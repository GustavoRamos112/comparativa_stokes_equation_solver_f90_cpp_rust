#!/usr/bin/env python3
"""Compare output data directories from bump solver simulations.

Usage:
    python check_output.py <dir_a> <dir_b> [options]

Compares all xy and uv files between two data directories and reports
differences. Useful for verifying C++ output against Python/Fortran reference.
"""

import argparse
import os
import sys


def read_data_file(filepath: str) -> list[tuple[float, float]]:
  """Read a xy or uv file, returning a list of (float, float) tuples."""
  values: list[tuple[float, float]] = []
  with open(filepath) as f:
    for lineno, line in enumerate(f, 1):
      line: str = line.strip()
      if not line:
        continue
      parts: list[str] = line.split()
      if len(parts) != 2:
        raise ValueError(f"{filepath}:{lineno}: expected 2 values, got {len(parts)}")
      values.append((float(parts[0]), float(parts[1])))
  return values


def compare_data_dirs(
  dir_a: str,
  dir_b: str,
  subdirs: tuple[str, str] | None = None,
  rtol: float = 1e-5,
  atol: float = 1e-8,
  verbose: bool = False,
) -> dict:
  """Compare output data from two simulation directories.

  Args:
      dir_a: Path to first (reference) data directory.
      dir_b: Path to second data directory to compare.
      subdirs: Subdirectories to compare, e.g. ('xy', 'uv').
               Defaults to ('xy', 'uv').
      rtol: Relative tolerance for floating-point comparison.
      atol: Absolute tolerance for floating-point comparison.
      verbose: If True, print detailed per-line differences.

  Returns:
      Dict with keys:
          'total_lines': total lines compared,
          'differing_lines': lines with values outside tolerance,
          'max_abs_diff': maximum absolute difference found,
          'max_rel_diff': maximum relative difference found,
          'missing_in_a': files in dir_b but not dir_a,
          'missing_in_b': files in dir_a but not dir_b,
          'shape_errors': files with mismatched line counts,
          'file_errors': files that could not be read,
          'all_match': True if no differences found.
  """
  if subdirs is None:
    subdirs = ("xy", "uv")

  report: dict = {
    "total_lines": 0,
    "differing_lines": 0,
    "max_abs_diff": 0.0,
    "max_rel_diff": 0.0,
    "missing_in_a": [],
    "missing_in_b": [],
    "shape_errors": [],
    "file_errors": [],
    "all_match": True,
  }

  for subdir in subdirs:
    path_a: str = os.path.join(dir_a, subdir)
    path_b: str = os.path.join(dir_b, subdir)

    if not os.path.isdir(path_a):
      report["missing_in_b"].append(path_a)
      report["all_match"] = False
      continue
    if not os.path.isdir(path_b):
      report["missing_in_a"].append(path_b)
      report["all_match"] = False
      continue

    files_a: set[str] = {f for f in os.listdir(path_a) if f.endswith(".txt")}
    files_b: set[str] = {f for f in os.listdir(path_b) if f.endswith(".txt")}

    only_b: set[str] = files_b - files_a
    for f in sorted(only_b):
      report["missing_in_a"].append(os.path.join(path_a, f))
      report["all_match"] = False

    only_a: set[str] = files_a - files_b
    for f in sorted(only_a):
      report["missing_in_b"].append(os.path.join(path_b, f))
      report["all_match"] = False

    common: list[str] = sorted(files_a & files_b)

    for filename in common:
      file_a: str = os.path.join(path_a, filename)
      file_b: str = os.path.join(path_b, filename)

      try:
        vals_a: list[tuple[float, float]] = read_data_file(file_a)
      except (OSError, ValueError) as e:
        report["file_errors"].append(str(e))
        report["all_match"] = False
        continue

      try:
        vals_b: list[tuple[float, float]] = read_data_file(file_b)
      except (OSError, ValueError) as e:
        report["file_errors"].append(str(e))
        report["all_match"] = False
        continue

      if len(vals_a) != len(vals_b):
        report["shape_errors"].append(
          f"{subdir}/{filename}: {len(vals_a)} lines vs {len(vals_b)} lines"
        )
        report["all_match"] = False
        continue

      for i, (va, vb) in enumerate(zip(vals_a, vals_b)):
        report["total_lines"] += 1
        line_diff_abs = max(abs(va[0] - vb[0]), abs(va[1] - vb[1]))
        denom = max(abs(va[0]), abs(va[1]), abs(vb[0]), abs(vb[1]))
        line_diff_rel = line_diff_abs / denom if denom > atol else line_diff_abs

        if line_diff_abs > atol and line_diff_rel > rtol:
          report["differing_lines"] += 1
          report["all_match"] = False
          if verbose:
            print(
              f"  {subdir}/{filename}:{i + 1}: "
              f"({va[0]:.6e},{va[1]:.6e}) vs "
              f"({vb[0]:.6e},{vb[1]:.6e})  "
              f"abs_diff={line_diff_abs:.2e}  "
              f"rel_diff={line_diff_rel:.2e}"
            )

        report["max_abs_diff"] = max(report["max_abs_diff"], line_diff_abs)
        report["max_rel_diff"] = max(report["max_rel_diff"], line_diff_rel)

  return report


def print_report(
  report: dict, 
  label_a: str = "A", 
  label_b: str = "B"
) -> None:
  """Print a human-readable summary of comparison results."""
  if report["all_match"]:
    print(f"OK: All {report['total_lines']} lines match within tolerance.")
    return

  print(f"Differences found between {label_a} and {label_b}:")
  print()

  if report["missing_in_a"]:
    print(f"  Missing in {label_a}:")
    for f in report["missing_in_a"]:
      print(f"    {f}")
    print()

  if report["missing_in_b"]:
    print(f"  Missing in {label_b}:")
    for f in report["missing_in_b"]:
      print(f"    {f}")
    print()

  if report["shape_errors"]:
    print("  Shape mismatches (different line counts):")
    for e in report["shape_errors"]:
      print(f"    {e}")
    print()

  if report["file_errors"]:
    print("  File read errors:")
    for e in report["file_errors"]:
      print(f"    {e}")
    print()

  if report["differing_lines"]:
    print(f"  {report['differing_lines']} / {report['total_lines']} lines differ.")
    print(f"  Max absolute difference: {report['max_abs_diff']:.2e}")
    print(f"  Max relative difference: {report['max_rel_diff']:.2e}")
    print()

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description="Compare bump solver output data directories."
  )
  parser.add_argument("dir_a", help="Reference data directory")
  parser.add_argument("dir_b", help="Data directory to compare")
  parser.add_argument(
    "--tol",
    type=float,
    default=1e-8,
    help="Absolute tolerance (default: 1e-8)",
  )
  parser.add_argument(
    "--rtol",
    type=float,
    default=1e-5,
    help="Relative tolerance (default: 1e-5)",
  )
  parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="Print detailed per-line differences",
  )
  parser.add_argument(
    "-q",
    "--quiet",
    action="store_true",
    help="Only print summary (exit code indicates result)",
  )

  args: argparse.Namespace = parser.parse_args()

  report: dict = compare_data_dirs(
    args.dir_a,
    args.dir_b,
    rtol=args.rtol,
    atol=args.tol,
    verbose=args.verbose,
  )

  if not args.quiet:
    print_report(report, args.dir_a, args.dir_b)

  sys.exit(0 if report["all_match"] else 1)
