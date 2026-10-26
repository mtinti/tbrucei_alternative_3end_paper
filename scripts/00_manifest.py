#!/usr/bin/env python3
"""Write the run manifest: inputs with checksums, options and software versions.

`start` creates the manifest before the first step; `finish` adds the outputs,
their checksums and the completion time. Every number in docs/ can then be
traced to one recorded run.
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def file_record(path):
    rec = {"path": path, "bytes": os.path.getsize(path)}
    # large alignment files are recorded by size and modification time only
    if rec["bytes"] < 2 << 30:
        rec["sha256"] = sha256(path)
    rec["modified"] = datetime.fromtimestamp(os.path.getmtime(path), timezone.utc).isoformat()
    return rec


def versions():
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    for mod in ("pandas", "numpy", "scipy", "sklearn", "statsmodels", "pysam", "matplotlib"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            out[mod] = "not installed"
    try:
        r = subprocess.run(["Rscript", "-e", 'cat(R.version.string, "\\n"); for (p in c("DRIMSeq", "stageR", "DEXSeq", "edgeR")) '
                            'cat(p, tryCatch(as.character(packageVersion(p)), error = function(e) "not installed"), "\\n")'],
                           capture_output=True, text=True, timeout=120)
        for line in r.stdout.strip().splitlines():
            k, _, v = line.partition(" ")
            out[k if k != "R" else "R"] = v if k != "R" else line
    except Exception as e:
        out["R"] = f"not run ({e})"
    try:
        out["git_commit"] = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        out["git_dirty"] = bool(subprocess.run(["git", "status", "--porcelain", "--", "scripts", "run_analysis.sh"],
                                               capture_output=True, text=True).stdout.strip())
    except Exception:
        pass
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("action", choices=["start", "finish"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--inputs", nargs="*", default=[])
    ap.add_argument("--options", nargs="*", default=[], help="key=value")
    ap.add_argument("--outputs", nargs="*", default=[])
    args = ap.parse_args()
    now = datetime.now(timezone.utc).isoformat()
    if args.action == "start":
        m = {"started": now, "status": "running", "command": " ".join(sys.argv[1:]),
             "options": dict((o.split("=", 1) if "=" in o else (o, "")) for o in args.options), "versions": versions(),
             "inputs": [file_record(p) for p in args.inputs if os.path.isfile(p)]}
    else:
        m = json.load(open(args.out))
        m.update(finished=now, status="completed",
                 outputs=[file_record(p) for p in args.outputs if os.path.isfile(p)])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    json.dump(m, open(args.out, "w"), indent=1)
    print(f"manifest {args.action}: {args.out}")


if __name__ == "__main__":
    main()
