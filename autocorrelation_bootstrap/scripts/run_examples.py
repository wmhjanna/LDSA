#!/usr/bin/env python3
"""Run and validate every bundled BCN example and summary-level reproduction."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
from build_engine import build

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    measurements = []
    started = time.perf_counter()
    build()
    measurements.append(dict(step="compile", seconds=time.perf_counter()-started))
    def run(label, command):
        t = time.perf_counter()
        subprocess.run(list(map(str, command)), cwd=root, check=True)
        measurements.append(dict(step=label, seconds=time.perf_counter()-t))
    for setting, label in [("auto", "Auto"), ("iid", "iid"), ("3", "3"), ("6", "6"), ("12", "12")]:
        run(f"BCN_{label}_1999", [args.rscript, "scripts/01_theilsen_block_bootstrap.R",
            "data/BCN_LDSA_MetaFull.csv", out / f"BCN_{label}_1999.csv", 1999, setting])
    run("BCN_cross_599", [args.rscript, "scripts/01_theilsen_block_bootstrap.R",
        "data/BCN_cross_metric_monthly.csv", out / "BCN_cross_599.csv", 599, "auto"])
    run("REML_47", [args.rscript, "scripts/02_reml_synthesis.R", "reference_results/site_trends_47.csv", out / "REML"])
    run("BCN_paired_9999_all_settings", [sys.executable, "scripts/03_paired_trend_difference.py",
        "data/BCN_paired_monthly.csv", out / "paired", "--all-settings"])
    run("multiplicity_143", [sys.executable, "scripts/04_adjust_multiplicity.py",
        "reference_results/paired_differences_143_all_settings.csv", out / "multiplicity_143.csv", "--expected-family-size", 143])
    run("verification", [sys.executable, "scripts/verify_results.py", "--results", out])
    metadata = dict(python=platform.python_version(), system=platform.system(), machine=platform.machine(),
                    total_seconds=time.perf_counter()-started, steps=measurements)
    (out / "execution_times.json").write_text(json.dumps(metadata, indent=2))

if __name__ == "__main__":
    main()
