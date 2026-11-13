# Skill: RDKit MCP Server (pharmacy use)

## When to use this skill

Use the RDKit MCP Server (upstream: tandemai-inc/rdkit-mcp-server, MIT)
or the pharmacy-tailored `pharmacy-mcp/server.py` in this repo whenever the
user asks for **computed molecular properties** — never guess numbers a
model could calculate instead:

- "What is the molecular weight / logP / TPSA of <drug>?"
- "Does <drug> pass Lipinski's Rule of Five?"
- "Is this SMILES valid? What is its canonical form / formula?"
- "How similar are <drug A> and <drug B> structurally?"
- "What is the core scaffold of <drug>?"
- Converting between SMILES / SDF / PDB / 2D images, or generating 3D
  coordinates for a molecule.

Do NOT use it for doses, indications, interactions, contraindications,
side effects or prices — it knows none of those. Use a clinical reference
for those questions.

## Servers

| Server | Entry point | Transport | Tools |
|---|---|---|---|
| Upstream RDKit MCP Server | `run_server.py` | stdio / sse / streamable_http | Full RDKit 2025.3.1 surface (`python list_tools.py` to list) |
| Pharmacy MCP (this repo) | `pharmacy-mcp/server.py` | stdio | 7 pharmacy tools (below) |

Start the pharmacy server with stdio (default):

```bash
python pharmacy-mcp/server.py
```

Claude Desktop config (`claude_desktop_config.json`):

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

## Pharmacy tools (`pharmacy-mcp/server.py`)

| Tool | Arguments | Returns |
|---|---|---|
| `validate_smiles` | `smiles: str` | `{valid, canonical_smiles, formula, molecular_weight}` or `{valid: false, error}` |
| `drug_descriptors` | `smiles: str` | MW, exact MW, formula, Crippen logP, TPSA, HBD, HBA, rotatable bonds, heavy atoms, FractionCSP3 |
| `lipinski_check` | `smiles: str` | Per-rule pass/fail, `lipinski_violations` (0–4), `lipinski_pass`, Veber (TPSA ≤ 140, rotatable bonds ≤ 10) |
| `fingerprint_similarity` | `smiles1, smiles2, radius=2, n_bits=2048` | Tanimoto similarity 0–1 |
| `standardize_smiles` | `smiles: str` | canonical SMILES, standardized SMILES, InChI, InChIKey |
| `murcko_scaffold` | `smiles: str` | Murcko scaffold SMILES |
| `drug_profile_by_name` | `drug_name: str` | descriptors + Lipinski for aspirin, caffeine, metformin, paracetamol, ibuprofen (sample table only) |

## Upstream tools worth knowing

The full server exposes the same RDKit functions directly, e.g.
`MolWt`, `ExactMolWt`, `CalcTPSA`, `CalcCrippenDescriptors` (logP),
`CalcNumHBA`, `CalcNumHBD`, `CalcNumRotatableBonds`, `CalcMolFormula`,
`TanimotoSimilarity`, `MurckoScaffoldSmilesFromSmiles`, `MolToImage`,
`EmbedMolecule`, plus format converters (`smiles_to_mol`,
`mol_to_smiles`, `mol_to_sdf`, `sdf_to_mol`, `mol_to_pdb`) and
`batch_map` for running one tool over a list of inputs. List them all with:

```bash
python list_tools.py
```

## Example calls

User: "Does aspirin pass Lipinski, and what are its logP and TPSA?"

1. `drug_descriptors(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")`
   → MW 180.16, formula C9H8O4, logP ≈ 1.31, TPSA ≈ 63.6
2. `lipinski_check(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")`
   → 0 violations, pass.

User: "How similar are aspirin and ibuprofen?"

1. `fingerprint_similarity(smiles1="CC(=O)OC1=CC=CC=C1C(=O)O", smiles2="CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")`
   → Tanimoto ≈ 0.2–0.3 (dissimilar — different scaffolds).

User: "What is the scaffold of metformin?"

1. `murcko_scaffold(smiles="CN(C)C(=N)NC(=N)N")` → scaffold of the biguanide core.

User: "Give me the profile of metformin."

1. `drug_profile_by_name(drug_name="metformin")`
   → descriptors + Lipinski from the sample table. **Always repeat the
   disclaimer to the user** (see below).

## Iranian-market rule (mandatory)

The `SAMPLE_DRUGS` table in `pharmacy-mcp/server.py` is a **small sample,
labelled sample/unverified**. It contains standard SMILES for a few
well-known generics and **international** originator brand names as
examples only. It asserts NOTHING about Iranian local brands, Iranian
availability, doses or prices.

- NEVER present a brand, dose, price or availability from that table as an
  Iranian-market fact.
- For any Iranian availability / dose / price question, verify against the
  user's own sources: the **Payamdaroosaz website** and/or the **IFDA
  yearbook**. If you cannot verify, say so explicitly.
- The sample table covers only aspirin, caffeine, metformin, paracetamol
  and ibuprofen. For any other drug, ask the user for its SMILES (or look
  it up in a trusted structure source) and use `drug_descriptors` /
  `lipinski_check` with that SMILES.

## Limits

- Computed logP (Crippen) and TPSA are **estimates**, not experimental values.
- Lipinski pass is a drug-likeness screen, not proof of oral
  bioavailability or clinical suitability.
- Wrong SMILES in → wrong numbers out. If a tool returns
  `valid: false` / an error, tell the user the SMILES is invalid and ask
  for a corrected one; do not invent a replacement SMILES.
- Tanimoto > 0.7 is often treated as "structurally similar", < 0.3 as
  dissimilar — a rule of thumb, not a clinical claim.
- This server gives no dose, interaction, indication or price
  information. Never use it to make a clinical decision.
