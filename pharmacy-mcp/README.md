# pharmacy-mcp — Pharmacy MCP Server for Iran

A small, standalone MCP server (Python, **stdio** transport) built on the
RDKit engine of the upstream
[RDKit MCP Server](https://github.com/tandemai-inc/rdkit-mcp-server)
(MIT, Copyright (c) 2025 TandemAI). It exposes 7 pharmacy-oriented tools
for AI agents: drug descriptors, Lipinski Rule-of-Five check, SMILES
validation/standardization, fingerprint similarity, Murcko scaffold, and a
small sample drug lookup.

Upstream license: MIT — see `../LICENSE`. Original attribution is kept;
nothing from the upstream repo was removed.

## Tools

| Tool | What it does |
|---|---|
| `validate_smiles` | Validate a SMILES string; returns canonical SMILES, formula, MW |
| `drug_descriptors` | MW, exact MW, formula, Crippen logP, TPSA, HBD/HBA, rotatable bonds, FractionCSP3 |
| `lipinski_check` | Lipinski Rule of Five (+ Veber) pass/fail per rule, violation count |
| `fingerprint_similarity` | Morgan-fingerprint Tanimoto similarity (0–1) between two SMILES |
| `standardize_smiles` | Canonical SMILES + InChI/InChIKey |
| `murcko_scaffold` | Murcko scaffold SMILES of a molecule |
| `drug_profile_by_name` | Descriptors + Lipinski for a drug in the **sample** table (aspirin, caffeine, metformin, paracetamol, ibuprofen) |

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Or reuse the upstream repo's environment (from the repo root):

```bash
pip install .
```

Requires Python >= 3.10 and RDKit 2025.3.1 (see `requirements.txt`).

> **Note on the MCP SDK version:** pin `mcp>=1.23.0,<2.0.0`. The upstream
> repo (and this server) import `FastMCP` from `mcp.server.fastmcp`; in
> MCP SDK 2.x that import was removed (`FastMCP` was renamed to
> `MCPServer`), so an unpinned `pip install .` currently pulls mcp 2.x and
> breaks both servers at import time. Pinning to 1.x is the fix until the
> upstream code is migrated.

## Run

stdio transport (default — for MCP clients such as Claude Desktop, or any
MCP-compatible agent):

```bash
python server.py
```

Example MCP client config (Claude Desktop `claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "pharmacy-mcp": {
      "command": "/absolute/path/to/.venv/bin/python",
      "args": ["/absolute/path/to/pharmacy-mcp/server.py"]
    }
  }
}
```

## Smoke test

```bash
python test_smoke.py
```

The smoke test imports the server module, calls the tool functions
directly (no MCP client needed) and checks descriptors for aspirin and
caffeine, the Lipinski check, a similarity call and an invalid-SMILES call.
It prints PASS/FAIL per check and exits non-zero on failure.

## Iranian market data — read this

The `SAMPLE_DRUGS` table in `server.py` is a **small sample, clearly
labelled sample/unverified**. It maps a few well-known generics to their
standard SMILES and to **international** originator brand names as examples
only. It asserts **nothing** about Iranian local brands, Iranian
availability, doses or prices.

Agents using this server MUST verify Iranian availability, dose and price
against the user's own sources — the **Payamdaroosaz website** and/or the
**IFDA yearbook** — and must never quote the sample table as Iranian market
fact. The server itself gives no dose, interaction, indication or price
information, and must not be used for clinical decisions.

## Limits

- Computed descriptors (logP, TPSA) are RDKit estimates, not experimental values.
- Lipinski pass is a drug-likeness screen, not a guarantee of oral
  bioavailability or clinical suitability.
- SMILES in → numbers out: a wrong SMILES gives wrong numbers.
- This server does **not** know doses, interactions, indications,
  contraindications or prices. Use a clinical reference for those.
