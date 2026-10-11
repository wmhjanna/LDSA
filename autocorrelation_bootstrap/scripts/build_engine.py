#!/usr/bin/env python3
"""Build the weighted Theil-Sen engine using a local C++ compiler."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys

def library_path():
    name = "weighted_sen.dylib" if sys.platform == "darwin" else "weighted_sen.dll" if os.name == "nt" else "weighted_sen.so"
    return Path(__file__).resolve().parent / name

def build(compiler=None):
    root = Path(__file__).resolve().parent
    compiler = shlex.split(compiler or os.environ.get("CXX", "c++"))
    flags = ["-dynamiclib"] if sys.platform == "darwin" else ["-shared", "-fPIC"] if os.name != "nt" else ["-shared"]
    output = library_path()
    subprocess.run(compiler + ["-O3", "-std=c++11"] + flags + [str(root / "weighted_sen.cpp"), "-o", str(output)], check=True)
    print(output)
    return output

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", help="Compiler command, default CXX or c++")
    build(parser.parse_args().compiler)
