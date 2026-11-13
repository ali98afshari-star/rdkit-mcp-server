"""Pharmacy MCP server for Iranian pharmacy needs.

Built on the RDKit engine of the upstream RDKit MCP Server
(tandemai-inc/rdkit-mcp-server, MIT, Copyright (c) 2025 TandemAI).
This file is a NEW, standalone FastMCP (stdio) server: it imports RDKit
directly and exposes a small set of pharmacy-oriented tools.

Tools:
  1. validate_smiles        - validate a SMILES string, return canonical SMILES + formula + MW
  2. drug_descriptors       - MW, logP (Crippen), TPSA, HBA, HBD, rotatable bonds, FractionCSP3
  3. lipinski_check         - Lipinski Rule of Five (+ Veber) pass/fail per rule
  4. fingerprint_similarity - Morgan-fingerprint Tanimoto similarity between two SMILES
  5. standardize_smiles     - canonical SMILES (and InChI when available)
  6. murcko_scaffold        - Murcko scaffold SMILES of a molecule
  7. drug_profile_by_name   - descriptors + Lipinski for a drug from the SAMPLE table

IMPORTANT - Iranian market data:
  The SAMPLE_DRUGS table below is a SMALL SAMPLE ONLY, clearly labelled
  sample/unverified. It maps a few well-known generic drugs to their
  standard SMILES and to INTERNATIONAL originator brand names as examples.
  It asserts NOTHING about Iranian local brands, Iranian availability,
  doses or prices. Agents MUST verify Iranian availability, dose and
  price against the user's own sources (Payamdaroosaz website /
  IFDA yearbook) and never quote this sample as Iranian market fact.

Run (stdio transport, the default for MCP clients):
    python server.py
"""

from __future__ import annotations

import json

from mcp.server.fastmcp import FastMCP
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

mcp = FastMCP("Pharmacy-MCP")

# ---------------------------------------------------------------------------
# SAMPLE table (UNVERIFIED for the Iranian market - see module docstring).
# SMILES strings are standard public chemistry data for these generics.
# "brand_example" is an international originator brand name, given only as
# an example of what a brand name looks like - it is NOT an Iranian brand
# claim and NOT a statement of Iranian availability.
# ---------------------------------------------------------------------------
SAMPLE_DRUGS = {
    "aspirin": {
        "generic": "aspirin",
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "brand_example": "Bayer Aspirin (international originator brand, SAMPLE)",
    },
    "caffeine": {
        "generic": "caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "brand_example": "No single originator brand (SAMPLE)",
    },
    "metformin": {
        "generic": "metformin",
        "smiles": "CN(C)C(=N)NC(=N)N",
        "brand_example": "Glucophage (international originator brand, SAMPLE)",
    },
    "paracetamol": {
        "generic": "paracetamol (acetaminophen)",
        "smiles": "CC(=O)NC1=CC=C(C=C1)O",
        "brand_example": "Tylenol / Panadol (international brands, SAMPLE)",
    },
    "ibuprofen": {
        "generic": "ibuprofen",
        "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "brand_example": "Brufen / Advil (international brands, SAMPLE)",
    },
}

SAMPLE_DISCLAIMER = (
    "SAMPLE DATA - UNVERIFIED for the Iranian market. "
    "Brand names above are international examples only. "
    "Verify Iranian availability, dose and price on the Payamdaroosaz "
    "website and/or the IFDA yearbook; never quote this sample as fact."
)


def _mol_or_error(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles!r}")
    return mol


def _descriptors(mol) -> dict:
    return {
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "exact_molecular_weight": round(Descriptors.ExactMolWt(mol), 4),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "logP_crippen": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "h_bond_donors": Lipinski.NumHDonors(mol),
        "h_bond_acceptors": Lipinski.NumHAcceptors(mol),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "fraction_csp3": round(Lipinski.FractionCSP3(mol), 2),
    }


@mcp.tool()
def validate_smiles(smiles: str) -> str:
    """Validate a SMILES string. Returns canonical SMILES, molecular formula
    and molecular weight, or an error if the SMILES is invalid."""
    try:
        mol = _mol_or_error(smiles)
    except ValueError as exc:
        return json.dumps({"valid": False, "error": str(exc)}, ensure_ascii=False)
    return json.dumps(
        {
            "valid": True,
            "input_smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
        },
        ensure_ascii=False,
    )


@mcp.tool()
def drug_descriptors(smiles: str) -> str:
    """Calculate drug-like descriptors for a molecule given as SMILES:
    MW, exact MW, formula, Crippen logP, TPSA, H-bond donors/acceptors,
    rotatable bonds, heavy atoms, FractionCSP3. Returns JSON."""
    try:
        mol = _mol_or_error(smiles)
    except ValueError as exc:
        return json.dumps({"error": str(exc)}, ensure_ascii=False)
    return json.dumps(_descriptors(mol), ensure_ascii=False)


@mcp.tool()
def lipinski_check(smiles: str) -> str:
    """Check Lipinski's Rule of Five (plus Veber rules) for a molecule given
    as SMILES. Returns each rule, its value, pass/fail, and the total number
    of Lipinski violations (0-4). Returns JSON."""
    try:
        mol = _mol_or_error(smiles)
    except ValueError as exc:
        return json.dumps({"error": str(exc)}, ensure_ascii=False)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    rules = [
        {"rule": "MW <= 500 Da", "value": round(mw, 2), "pass": mw <= 500},
        {"rule": "logP <= 5", "value": round(logp, 2), "pass": logp <= 5},
        {"rule": "H-bond donors <= 5", "value": hbd, "pass": hbd <= 5},
        {"rule": "H-bond acceptors <= 10", "value": hba, "pass": hba <= 10},
    ]
    violations = sum(1 for r in rules if not r["pass"])
    return json.dumps(
        {
            "rules": rules,
            "lipinski_violations": violations,
            "lipinski_pass": violations <= 1,
            "veber": {
                "tpsa <= 140 A^2": {"value": round(tpsa, 2), "pass": tpsa <= 140},
                "rotatable_bonds <= 10": {"value": rotb, "pass": rotb <= 10},
            },
            "note": "Lipinski pass is a drug-likeness screen, not a guarantee "
            "of oral bioavailability or clinical suitability.",
        },
        ensure_ascii=False,
    )


@mcp.tool()
def fingerprint_similarity(smiles1: str, smiles2: str, radius: int = 2, n_bits: int = 2048) -> str:
    """Tanimoto similarity (0-1) between two molecules given as SMILES,
    using Morgan fingerprints. ~1 = near-identical, >0.7 often considered
    structurally similar, <0.3 dissimilar. Returns JSON."""
    try:
        mol1 = _mol_or_error(smiles1)
        mol2 = _mol_or_error(smiles2)
    except ValueError as exc:
        return json.dumps({"error": str(exc)}, ensure_ascii=False)
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=n_bits)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=n_bits)
    sim = AllChem.DataStructs.TanimotoSimilarity(fp1, fp2)
    return json.dumps(
        {
            "smiles1": Chem.MolToSmiles(mol1),
            "smiles2": Chem.MolToSmiles(mol2),
            "tanimoto_similarity": round(sim, 4),
            "radius": radius,
            "n_bits": n_bits,
        },
        ensure_ascii=False,
    )


@mcp.tool()
def standardize_smiles(smiles: str) -> str:
    """Standardize a SMILES string: returns the canonical SMILES and the
    InChI/InChIKey when RDKit can generate them. Returns JSON."""
    try:
        mol = _mol_or_error(smiles)
    except ValueError as exc:
        return json.dumps({"error": str(exc)}, ensure_ascii=False)
    out = {
        "input_smiles": smiles,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
    }
    try:
        from rdkit.Chem import rdMolDescriptors as _rdmd  # noqa: F401
        from rdkit.Chem.MolStandardize import rdMolStandardize

        std = rdMolStandardize.StandardizeMol(mol)
        out["standardized_smiles"] = Chem.MolToSmiles(std)
    except Exception:
        out["standardized_smiles"] = out["canonical_smiles"]
    try:
        from rdkit.Chem import inchi

        out["inchi"] = inchi.MolToInchi(mol)
        out["inchi_key"] = inchi.MolToInchiKey(mol)
    except Exception:
        out["inchi"] = None
        out["inchi_key"] = None
    return json.dumps(out, ensure_ascii=False)


@mcp.tool()
def murcko_scaffold(smiles: str) -> str:
    """Extract the Murcko scaffold (core framework) of a molecule given as
    SMILES. Useful for grouping drugs by their core structure. Returns JSON."""
    try:
        mol = _mol_or_error(smiles)
    except ValueError as exc:
        return json.dumps({"error": str(exc)}, ensure_ascii=False)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return json.dumps(
        {
            "input_smiles": Chem.MolToSmiles(mol),
            "murcko_scaffold_smiles": Chem.MolToSmiles(scaffold),
        },
        ensure_ascii=False,
    )


@mcp.tool()
def drug_profile_by_name(drug_name: str) -> str:
    """Drug profile for a generic drug name from a SMALL SAMPLE table
    (aspirin, caffeine, metformin, paracetamol, ibuprofen): returns its
    SMILES, descriptors and Lipinski check. The table is SAMPLE/UNVERIFIED
    for the Iranian market - see the disclaimer in the result. For any drug
    not in the table, use validate_smiles/drug_descriptors with its SMILES
    instead. Returns JSON."""
    key = drug_name.strip().lower()
    entry = SAMPLE_DRUGS.get(key)
    if entry is None:
        return json.dumps(
            {
                "error": f"Drug {drug_name!r} is not in the sample table.",
                "sample_drugs": sorted(SAMPLE_DRUGS),
                "hint": "Use drug_descriptors / lipinski_check with the drug's SMILES instead.",
                "disclaimer": SAMPLE_DISCLAIMER,
            },
            ensure_ascii=False,
        )
    mol = Chem.MolFromSmiles(entry["smiles"])
    desc = _descriptors(mol)
    mw, logp = desc["molecular_weight"], desc["logP_crippen"]
    hbd, hba = desc["h_bond_donors"], desc["h_bond_acceptors"]
    violations = sum(
        [mw > 500, logp > 5, hbd > 5, hba > 10]
    )
    return json.dumps(
        {
            "generic": entry["generic"],
            "smiles": entry["smiles"],
            "brand_example": entry["brand_example"],
            "descriptors": desc,
            "lipinski_violations": violations,
            "lipinski_pass": violations <= 1,
            "disclaimer": SAMPLE_DISCLAIMER,
        },
        ensure_ascii=False,
    )


if __name__ == "__main__":
    mcp.run(transport="stdio")
