"""Collect completed Phase-5 MCNP6.3 validation outputs.

Run ``run_mcnp_validation.cmd`` on the institutional computer. After copying
all ``mcnp\\sim_XXXX`` folders back, run this script on
the thesis computer with ``--collect``.  It never invokes MCNP and never
modifies an input deck.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path


CELL20 = re.compile(r"cell\s+20\s*\n\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s+([0-9.]+)", re.I)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_output(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="latin-1", errors="replace")
    low = text.lower()
    match = CELL20.search(text)
    nps = bool(re.search(r"run terminated when\s+100000000\s+particle histories were done", text, re.I))
    valid = bool(
        match and nps and "fatal error" not in low and "bad trouble" not in low
        and "you need a source subroutine" not in low and "passed all bin error check" in low
        and "tally        4 tfc bin did not pass" not in low and float(match.group(2)) < 0.10
    )
    return {
        "status_execucao_mcnp": "COMPLETED" if valid else "FAILED_REVIEW",
        "status_extracao": "COMPLETED" if valid else "FAILED_REVIEW",
        "arquivo_output": str(path), "runtime_s": "",
        "dose_efetiva_mSv_seg": "" if not match else format(float(match.group(1)), ".12g"),
        "erro_relativo": "" if not match else format(float(match.group(2)), ".12g"),
        "qa_tally": "PASS" if valid else "FAIL",
        "error_message": "" if valid else "Falha de término/tally/QA MCNP; consulte output.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-root", required=True, type=Path)
    parser.add_argument("--collect", action="store_true", help="Obrigatório: extrai os outputs já executados.")
    args = parser.parse_args()
    if not args.collect:
        raise SystemExit("Use --collect após copiar todos os outputs MCNP de volta.")
    source = args.study_root / "preparation_manifest.csv"
    rows = list(csv.DictReader(source.open(encoding="utf-8-sig", newline=""))) if source.is_file() else []
    if not rows:
        raise ValueError("Manifesto de preparação vazio ou ausente.")
    results = []
    for row in rows:
        sim_id = int(row["validation_sim_id"])
        inp = Path(row["mcnp_input"])
        if not inp.is_file() or sha256(inp).lower() != row["mcnp_input_sha256"].lower():
            raise ValueError(f"Input MCNP ausente ou com hash divergente: {inp}")
        output = inp.parent / f"output_{sim_id}.o"
        if not output.is_file():
            raise FileNotFoundError(output)
        results.append({
            "validation_candidate_id": row["validation_candidate_id"],
            "validation_sim_id": sim_id, "mcnp_input_sha256": row["mcnp_input_sha256"],
            **parse_output(output),
        })
    destination = args.study_root / "mcnp_run_manifest.csv"
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader(); writer.writerows(results)
    if any(row["qa_tally"] != "PASS" for row in results):
        raise RuntimeError(f"Coleta concluída com falha; consulte {destination}")
    print(f"UTC={datetime.now(timezone.utc).isoformat()}")
    print(destination)


if __name__ == "__main__":
    main()
