"""Smoke test for pharmacy-mcp/server.py.

Runs the tool functions directly (no MCP client needed) and checks:
  1. aspirin descriptors (MW ~180.16, formula C9H8O4)
  2. caffeine Lipinski check (0 violations)
  3. fingerprint similarity aspirin vs aspirin == 1.0
  4. invalid SMILES is rejected gracefully
  5. drug_profile_by_name('metformin') returns the sample table entry
     with the UNVERIFIED disclaimer
  6. standardize_smiles returns a canonical SMILES
  7. murcko_scaffold returns a non-empty scaffold

Usage:  python test_smoke.py
Exit code 0 = all checks PASS, 1 = any FAIL.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import server  # noqa: E402


def _call(tool, *args, **kwargs):
    """Call a @mcp.tool() function whether the decorator returned the plain
    function or a Tool wrapper exposing .fn."""
    fn = getattr(tool, "fn", tool)
    return fn(*args, **kwargs)


def _check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    return ok


def main() -> int:
    results = []

    # 1. Aspirin descriptors
    d = json.loads(_call(server.drug_descriptors, "CC(=O)OC1=CC=CC=C1C(=O)O"))
    results.append(_check(
        "aspirin descriptors",
        abs(d.get("molecular_weight", 0) - 180.16) < 0.05
        and d.get("formula") == "C9H8O4",
        f"MW={d.get('molecular_weight')} formula={d.get('formula')} "
        f"logP={d.get('logP_crippen')} TPSA={d.get('tpsa')}",
    ))

    # 2. Caffeine Lipinski check (0 violations expected)
    lip = json.loads(_call(server.lipinski_check, "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"))
    results.append(_check(
        "caffeine lipinski_check",
        lip.get("lipinski_violations") == 0 and lip.get("lipinski_pass") is True,
        f"violations={lip.get('lipinski_violations')}",
    ))

    # 3. Self-similarity == 1.0
    sim = json.loads(_call(
        server.fingerprint_similarity,
        "CC(=O)OC1=CC=CC=C1C(=O)O", "CC(=O)OC1=CC=CC=C1C(=O)O",
    ))
    results.append(_check(
        "fingerprint_similarity self == 1.0",
        abs(sim.get("tanimoto_similarity", 0) - 1.0) < 1e-6,
        f"sim={sim.get('tanimoto_similarity')}",
    ))

    # 4. Invalid SMILES rejected gracefully
    bad = json.loads(_call(server.validate_smiles, "not_a_smiles!!!"))
    results.append(_check(
        "invalid SMILES rejected",
        bad.get("valid") is False,
        f"valid={bad.get('valid')}",
    ))

    # 5. Sample drug profile (metformin) + disclaimer present
    prof = json.loads(_call(server.drug_profile_by_name, "metformin"))
    results.append(_check(
        "drug_profile_by_name metformin",
        prof.get("generic") == "metformin"
        and "UNVERIFIED" in prof.get("disclaimer", ""),
        f"MW={prof.get('descriptors', {}).get('molecular_weight')}",
    ))

    # 6. Standardize SMILES
    std = json.loads(_call(server.standardize_smiles, "CC(=O)OC1=CC=CC=C1C(=O)O"))
    results.append(_check(
        "standardize_smiles",
        bool(std.get("canonical_smiles")),
        f"canonical={std.get('canonical_smiles')}",
    ))

    # 7. Murcko scaffold
    sc = json.loads(_call(server.murcko_scaffold, "CC(=O)OC1=CC=CC=C1C(=O)O"))
    results.append(_check(
        "murcko_scaffold",
        bool(sc.get("murcko_scaffold_smiles")),
        f"scaffold={sc.get('murcko_scaffold_smiles')}",
    ))

    print(f"\n{sum(results)}/{len(results)} checks passed.")
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
