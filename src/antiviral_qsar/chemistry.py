from __future__ import annotations

import math
import re
from decimal import Decimal, InvalidOperation
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from .schema import Measurement

RDLogger.DisableLog("rdApp.error")


def structure(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("invalid_smiles")
    fragments = Chem.GetMolFrags(mol, asMols=True)
    organic = [m for m in fragments if any(a.GetAtomicNum() == 6 for a in m.GetAtoms())]
    if len(organic) > 1:
        raise ValueError("multiple_organic_fragments")
    if not organic:
        raise ValueError("not_small_organic_molecule")
    parent = rdMolStandardize.Cleanup(organic[0])
    if parent.GetNumHeavyAtoms() > 150:
        raise ValueError("outside_small_molecule_scope")
    canonical = Chem.MolToSmiles(parent, isomericSmiles=True, canonical=True)
    return canonical, Chem.MolToInchiKey(parent), len(fragments) > 1


def numeric_value(value: str, relation: str = "="):
    raw = str(value).strip().replace("−", "-").replace("≤", "<=").replace("≥", ">=")
    match = re.fullmatch(r"\s*(<=|>=|<|>|=|~|≈)?\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*", raw)
    if not match:
        raise ValueError("non_scalar_value")
    prefix, number = match.groups()
    relation = relation.strip().replace("≤", "<=").replace("≥", ">=") or "="
    if prefix:
        if relation not in ("=", prefix):
            raise ValueError("conflicting_relation")
        relation = prefix
    if relation not in ("=", ">", "<", ">=", "<="):
        raise ValueError("unsupported_relation")
    try:
        v = Decimal(number)
    except InvalidOperation:
        raise ValueError("invalid_number") from None
    if not v.is_finite() or v <= 0:
        raise ValueError("nonpositive_or_nonfinite_value")
    return v, relation


def normalize(value: str, units: str, relation="="):
    v, rel = numeric_value(value, relation)
    unit = units.strip().replace("µ", "u").replace("μ", "u").replace(" ", "")
    factors = {"M": "1e9", "mM": "1e6", "uM": "1e3", "nM": "1", "pM": "1e-3", "fM":"1e-6",
               "mol/L": "1e9", "mmol/L": "1e6", "umol/L": "1e3", "nmol/L": "1"}
    if unit not in factors:
        raise ValueError("unsupported_units")
    nm = float(v * Decimal(factors[unit]))
    if not math.isfinite(nm) or nm <= 0:
        raise ValueError("normalized_value_out_of_range")
    return nm, rel, 9 - math.log10(nm) if rel == "=" else None


def validate_measurement(m: Measurement, *, llm_gate=False):
    reasons = list(m.review_reasons)
    if m.smiles_original:
        try:
            m.smiles, m.inchikey, salt = structure(m.smiles_original)
            if salt:
                m.warnings.append("counterions_removed_original_preserved")
        except ValueError as e:
            reasons.append(str(e))
    else:
        reasons.append("unresolved_structure")
    try:
        m.value_nm, m.relation, m.p_activity = normalize(m.value_original, m.units_original, m.relation)
    except ValueError as e:
        reasons.append(str(e))
    if m.endpoint not in ("IC50", "EC50", "CC50"):
        reasons.append("unsupported_endpoint")
    if not m.cell_line:
        reasons.append("missing_cell_line")
    if not m.cell_context_supported:
        reasons.append("cellular_context_not_supported")
    if m.endpoint != "CC50" and not m.virus:
        reasons.append("missing_virus")
    if not m.evidence_grounded:
        reasons.append("evidence_not_grounded")
    if m.extraction_method == "llm" and not llm_gate:
        reasons.append("extractor_not_validated")
    if m.extraction_method == "structured" and m.context_method == "llm":
        reasons.append("assay_classifier_requires_review")
    m.review_reasons = sorted(set(reasons))
    m.warnings = sorted(set(m.warnings))
    return m
