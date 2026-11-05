from __future__ import annotations

from typing import Literal
from pydantic import BaseModel, Field, ConfigDict


class Virus(BaseModel):
    query: str
    name: str
    tax_id: str
    synonyms: list[str] = Field(default_factory=list)
    scope: str = "selected taxon and its subtypes; individual strains retained"


class AssayContext(BaseModel):
    model_config = ConfigDict(extra="ignore")
    kind: Literal["cellular_antiviral", "cytotoxicity", "biochemical", "other", "unclear"] = "unclear"
    cell_line: str | None = None
    virus: str | None = None
    strain: str | None = None
    duration: str | None = None
    method: str | None = None
    evidence: str = ""
    supported: bool = False
    classified_by: str = "rules"


class ExtractedMeasurement(BaseModel):
    model_config = ConfigDict(extra="forbid")
    compound: str
    smiles: str | None = None
    endpoint: Literal["IC50", "EC50", "CC50"]
    value: str
    units: str
    relation: str = "="
    virus: str | None = None
    cell_line: str | None = None
    strain: str | None = None
    duration: str | None = None
    method: str | None = None
    experiment_label: str | None = None
    evidence: str
    context_evidence: str = ""


class Extraction(BaseModel):
    model_config = ConfigDict(extra="forbid")
    measurements: list[ExtractedMeasurement] = Field(default_factory=list)


class Measurement(BaseModel):
    observation_id: str
    source: str
    source_record_id: str
    source_url: str
    source_snapshot: str = ""
    source_assay_id: str = ""
    origin_assay_id: str = ""
    document_id: str = ""
    doi: str = ""
    pmids: list[str] = Field(default_factory=list)
    compound: str = ""
    compound_id: str = ""
    smiles_original: str = ""
    smiles: str = ""
    inchikey: str = ""
    virus: str = ""
    virus_tax_id: str = ""
    strain: str = ""
    cell_line: str = ""
    endpoint: str = ""
    value_original: str = ""
    units_original: str = ""
    relation: str = "="
    value_nm: float | None = None
    p_activity: float | None = None
    duration: str = ""
    method: str = ""
    assay_description: str = ""
    evidence: str = ""
    evidence_locator: str = ""
    context_evidence: str = ""
    extraction_method: str = "structured"
    context_method: str = "rules"
    extractor_signature: str = ""
    evidence_grounded: bool = False
    cell_context_supported: bool = False
    review_reasons: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    provenance: list[dict] = Field(default_factory=list)
    duplicate_of: str = ""
    possible_duplicate_group: str = ""
    literature_role: str = ""
    supplement_of: list[str] = Field(default_factory=list)

    @property
    def status(self):
        return "review" if self.review_reasons else "accepted"
