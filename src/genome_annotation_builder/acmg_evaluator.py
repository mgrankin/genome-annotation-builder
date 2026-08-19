"""Conservative, disease-aware ACMG secondary-findings candidate evaluation."""

from __future__ import annotations

import hashlib
import re
from collections import defaultdict
from datetime import UTC, datetime
from typing import Any

from .constants import ACMG_SF_GENE_SET_SHA256, ACMG_SF_SNAPSHOT, ACMG_SF_VERSION
from .util import canonical_json_bytes

PATHOGENIC_SIGNIFICANCE = {"Pathogenic", "Likely_pathogenic"}
BENIGN_SIGNIFICANCE = {"Benign", "Likely_benign"}
REPORTABLE_CURATION = {"Yes", "See Guidance"}
TRUNCATING_TERMS = {
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_gained",
}
MECHANISM_REVIEW_GENES = {"PCSK9", "RET", "TGFBR1"}


def _policy_key(row: dict[str, Any]) -> str:
    payload = {
        "gene_symbol": row["gene_symbol"],
        "disease": row["disease"],
        "inheritance": row["inheritance"],
        "disorder_mim_ids": row["disorder_mim_ids"],
    }
    return hashlib.sha256(canonical_json_bytes(payload)).hexdigest()


def _condition_ids(values: list[str]) -> set[str]:
    output: set[str] = set()
    for value in values:
        for prefix, identifier in re.findall(r"(?i)(MONDO(?::MONDO)?|OMIM|MIM)\s*:\s*(\d+)", value):
            canonical_prefix = "MONDO" if prefix.upper().startswith("MONDO") else "MIM"
            output.add(f"{canonical_prefix}:{identifier}")
    return output


def _gene_tokens(values: list[str]) -> set[str]:
    return {value.split(":", maxsplit=1)[0].strip().upper() for value in values if value.strip()}


def _disease_tokens(value: str) -> set[str]:
    stop = {
        "a",
        "an",
        "and",
        "autosomal",
        "disease",
        "disorder",
        "hereditary",
        "of",
        "related",
        "syndrome",
        "the",
        "type",
    }
    return {
        token
        for token in re.findall(r"[a-z0-9]+", value.lower())
        if token not in stop and len(token) > 1
    }


def _transcript_terms(transcripts: list[dict[str, Any]]) -> set[str]:
    return {
        str(term)
        for transcript in transcripts
        for term in transcript.get("consequence_terms") or []
    }


def _hgvsp_values(transcripts: list[dict[str, Any]]) -> set[str]:
    return {str(transcript["hgvsp"]) for transcript in transcripts if transcript.get("hgvsp")}


def _trdn_transcript_in_scope(transcripts: list[dict[str, Any]]) -> bool:
    for transcript in transcripts:
        refseq = str(transcript.get("refseq_transcript_id") or "")
        transcript_id = str(transcript.get("transcript_id") or "")
        if not (refseq.startswith("NM_001256021") or transcript_id.startswith("NM_001256021")):
            continue
        exon = str(transcript.get("exon") or "").split("/", maxsplit=1)[0]
        if exon.isdigit() and 1 <= int(exon) <= 8:
            return True
    return False


def _snapshot_datetime(snapshot: str) -> datetime:
    try:
        return datetime.strptime(snapshot, "%Y%m%d").replace(tzinfo=UTC)
    except ValueError as exc:
        raise ValueError(f"invalid ACMG/ClinGen snapshot date: {snapshot!r}") from exc


class AcmgSecondaryFindingsEvaluator:
    """Accumulate a small set of policy candidates and resolve genotype rules."""

    def __init__(self, snapshot: dict[str, Any], min_review_stars: int = 2) -> None:
        self.version = str(snapshot["version"])
        self.snapshot = str(snapshot["snapshot"])
        if self.version != ACMG_SF_VERSION or self.snapshot != ACMG_SF_SNAPSHOT:
            raise ValueError("ACMG evaluator received an unpinned policy snapshot")
        self.min_review_stars = min_review_stars
        self.policy_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in snapshot["policy_gene_disease_pairs"]:
            policy = dict(row)
            policy["policy_key"] = _policy_key(policy)
            self.policy_by_gene[str(policy["gene_symbol"]).upper()].append(policy)
        self.genes = frozenset(self.policy_by_gene)
        self.guidance = {
            str(row["gene_symbol"]).upper(): row.get("comment")
            for row in snapshot["clingen_reporting_guidance"]
        }
        self.target_ids, self.curation_status = self._build_disease_targets(
            snapshot["clingen_gene_disease_reporting"]
        )
        self._decisions: list[dict[str, Any]] = []

    def _build_disease_targets(
        self, disease_rows: list[dict[str, Any]]
    ) -> tuple[dict[str, set[str]], dict[str, list[str]]]:
        target_ids: dict[str, set[str]] = {}
        status: dict[str, list[str]] = {}
        by_group: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
        for row in disease_rows:
            by_group[(str(row["gene_symbol"]).upper(), str(row["mode_of_inheritance"]))].append(row)

        for gene, policies in self.policy_by_gene.items():
            grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
            for policy in policies:
                grouped[str(policy["inheritance"])].append(policy)
            for inheritance, same_inheritance in grouped.items():
                curation = by_group.get((gene, inheritance), [])
                reportable = [
                    row
                    for row in curation
                    if row["reportable_as_secondary_finding"] in REPORTABLE_CURATION
                ]
                assignments: dict[str, list[dict[str, Any]]] = defaultdict(list)
                if len(same_inheritance) == 1:
                    assignments[same_inheritance[0]["policy_key"]] = reportable
                else:
                    for disease in reportable:
                        disease_terms = _disease_tokens(str(disease["disease"]))
                        scores = [
                            len(disease_terms & _disease_tokens(str(policy["disease"])))
                            for policy in same_inheritance
                        ]
                        best = max(scores, default=0)
                        if best:
                            for policy, score in zip(same_inheritance, scores, strict=True):
                                if score == best:
                                    assignments[policy["policy_key"]].append(disease)
                for policy in same_inheritance:
                    key = str(policy["policy_key"])
                    targets = {str(value) for value in policy["disorder_mim_ids"]}
                    targets.update(
                        str(row["mondo_id"])
                        for row in assignments.get(key, [])
                        if row.get("mondo_id")
                    )
                    target_ids[key] = targets
                    status[key] = sorted(
                        {
                            str(row["reportable_as_secondary_finding"])
                            for row in assignments.get(key, [])
                        }
                    )
        return target_ids, status

    def observe(
        self,
        *,
        variant_id: str,
        gene_symbols: set[str],
        gt: tuple[int, ...],
        phased: bool,
        call_confidence: str,
        transcripts: list[dict[str, Any]],
        clinvar_rows: list[dict[str, Any]],
    ) -> None:
        """Add exact-allele evidence; unrelated/benign evidence is retained elsewhere."""

        for gene in sorted({symbol.upper() for symbol in gene_symbols} & self.genes):
            gene_rows = [
                row
                for row in clinvar_rows
                if gene in _gene_tokens([str(value) for value in row.get("gene_info", [])])
            ]
            if not gene_rows:
                continue
            for policy in self.policy_by_gene[gene]:
                targets = self.target_ids[str(policy["policy_key"])]
                matched: list[tuple[dict[str, Any], set[str]]] = []
                for evidence in gene_rows:
                    identifiers = _condition_ids(
                        [str(value) for value in evidence.get("condition_identifiers", [])]
                    )
                    overlap = identifiers & targets
                    if overlap:
                        matched.append((evidence, overlap))
                if not matched:
                    continue
                significance = {str(row.get("aggregate_significance")) for row, _overlap in matched}
                has_pathogenic = bool(significance & PATHOGENIC_SIGNIFICANCE)
                has_benign = bool(significance & BENIGN_SIGNIFICANCE)
                has_conflicts = (has_pathogenic and has_benign) or any(
                    bool(row.get("has_conflicts")) for row, _overlap in matched
                )
                if not has_pathogenic and not has_conflicts:
                    continue
                pathogenic_rows = [
                    row
                    for row, _overlap in matched
                    if row.get("aggregate_significance") in PATHOGENIC_SIGNIFICANCE
                ]
                stars = max((int(row.get("review_stars", 0)) for row in pathogenic_rows), default=0)
                aggregate = (
                    "Pathogenic"
                    if any(
                        row.get("aggregate_significance") == "Pathogenic" for row in pathogenic_rows
                    )
                    else "Likely_pathogenic"
                    if pathogenic_rows
                    else "Conflicting_interpretations"
                )
                clinvar_snapshots = {str(row["snapshot"]) for row, _overlap in matched}
                if len(clinvar_snapshots) != 1:
                    raise ValueError(
                        "ACMG evaluation requires one immutable ClinVar snapshot per build"
                    )
                decision = {
                    "variant_id": variant_id,
                    "gene_symbol": gene,
                    "policy_key": policy["policy_key"],
                    "disease": policy["disease"],
                    "disease_identifiers": sorted(targets),
                    "inheritance": policy["inheritance"],
                    "phenotype_category": policy["phenotype_category"],
                    "variants_to_report": policy["variants_to_report"],
                    "decision_status": "PENDING_POLICY_EVALUATION",
                    "is_actionable": False,
                    "clinical_grade": False,
                    "confirmation_required": True,
                    "analytical_confidence": call_confidence,
                    "genotype": list(gt),
                    "phased": phased,
                    "clinvar_significance": aggregate,
                    "clinvar_snapshot": next(iter(clinvar_snapshots)),
                    "clinvar_review_stars": stars,
                    "clinvar_has_conflicts": has_conflicts,
                    "clinvar_evidence_ids": sorted(
                        {str(row["evidence_id"]) for row, _overlap in matched}
                    ),
                    "clinvar_accessions": sorted(
                        {str(row["clinvar_accession"]) for row, _overlap in matched}
                    ),
                    "matched_condition_identifiers": sorted(
                        {identifier for _row, overlap in matched for identifier in overlap}
                    ),
                    "policy_version": self.version,
                    "policy_snapshot": self.snapshot,
                    "policy_gene_set_sha256": ACMG_SF_GENE_SET_SHA256,
                    "guidance": self.guidance.get(gene),
                    "reasons": [],
                    "_transcripts": transcripts,
                }
                self._apply_variant_rule(decision)
                self._decisions.append(decision)

    def _apply_variant_rule(self, decision: dict[str, Any]) -> None:
        if decision["clinvar_has_conflicts"]:
            decision["decision_status"] = "CLINVAR_CONFLICT"
            decision["reasons"].append("conflicting exact-allele ClinVar evidence")
            return
        if decision["clinvar_review_stars"] < self.min_review_stars:
            decision["decision_status"] = "INSUFFICIENT_CLINVAR_REVIEW"
            decision["reasons"].append(
                f"ClinVar review stars below configured minimum {self.min_review_stars}"
            )
            return
        if not self.curation_status.get(str(decision["policy_key"])):
            decision["decision_status"] = "DISEASE_REPORTABILITY_NOT_CONFIRMED"
            decision["reasons"].append("no reportable ClinGen disease row mapped to policy pair")
            return
        if decision["analytical_confidence"] != "high":
            decision["decision_status"] = "LOW_ANALYTICAL_CONFIDENCE"
            decision["reasons"].append("variant call does not meet configured high-confidence gate")
            return

        gene = str(decision["gene_symbol"])
        gt = tuple(int(value) for value in decision["genotype"])
        transcripts = decision["_transcripts"]
        terms = _transcript_terms(transcripts)
        if gene == "HFE":
            c282y = any(
                value.endswith(":p.Cys282Tyr") or value == "p.Cys282Tyr"
                for value in _hgvsp_values(transcripts)
            )
            if not c282y:
                decision["decision_status"] = "VARIANT_CLASS_NOT_REPORTABLE"
                decision["reasons"].append("ACMG SF reports only HFE p.Cys282Tyr")
            elif gt == (1, 1):
                self._mark_reportable(decision, "homozygous HFE p.Cys282Tyr")
            else:
                decision["decision_status"] = "GENOTYPE_NOT_REPORTABLE"
                decision["reasons"].append("HFE p.Cys282Tyr is not homozygous")
            return
        if gene == "APOB" and "missense_variant" not in terms:
            decision["decision_status"] = "VARIANT_CLASS_NOT_REPORTABLE"
            decision["reasons"].append(
                "APOB guidance requires a disease-associated missense variant"
            )
            return
        if gene == "TTN":
            if not terms & TRUNCATING_TERMS:
                decision["decision_status"] = "VARIANT_CLASS_NOT_REPORTABLE"
                decision["reasons"].append("TTN policy is restricted to truncating variants")
            else:
                decision["decision_status"] = "GENE_SPECIFIC_REVIEW_REQUIRED"
                decision["reasons"].append("TTN high-PSI exon applicability is not yet established")
            return
        if gene in MECHANISM_REVIEW_GENES:
            decision["decision_status"] = "GENE_SPECIFIC_REVIEW_REQUIRED"
            decision["reasons"].append(
                "gene-specific mechanism/domain applicability requires variant-level review"
            )
            return
        if gene == "TRDN" and not _trdn_transcript_in_scope(transcripts):
            decision["decision_status"] = "TRANSCRIPT_SCOPE_NOT_CONFIRMED"
            decision["reasons"].append("TRDN cardiac transcript exons 1-8 were not confirmed")
            return

        inheritance = str(decision["inheritance"])
        if inheritance == "AR":
            if gt == (1, 1):
                self._mark_reportable(decision, "homozygous qualifying allele")
            elif len(gt) == 2 and sum(gt) == 1:
                decision["decision_status"] = "PENDING_BIALLELIC_EVALUATION"
                decision["reasons"].append("a second qualifying allele and phase must be assessed")
            else:
                decision["decision_status"] = "GENOTYPE_NOT_REPORTABLE"
                decision["reasons"].append("genotype is not a supported biallelic state")
            return
        if gene == "ABCD1" and len(gt) == 2 and sum(gt) == 1:
            decision["decision_status"] = "PENDING_BIALLELIC_EVALUATION"
            decision["reasons"].append("a diploid heterozygote requires a second qualifying allele")
            return
        self._mark_reportable(decision, "policy, disease, evidence, and genotype gates passed")

    @staticmethod
    def _mark_reportable(decision: dict[str, Any], reason: str) -> None:
        decision["decision_status"] = "REPORTABLE_CANDIDATE"
        decision["is_actionable"] = True
        decision["reasons"].append(reason)

    def finalize(self) -> dict[str, Any]:
        pending: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for decision in self._decisions:
            if decision["decision_status"] == "PENDING_BIALLELIC_EVALUATION":
                pending[str(decision["policy_key"])].append(decision)
        for decisions in pending.values():
            distinct = {str(decision["variant_id"]) for decision in decisions}
            if len(distinct) >= 2:
                for decision in decisions:
                    decision["decision_status"] = "PHASE_REVIEW_REQUIRED"
                    decision["reasons"].append(
                        "multiple qualifying heterozygous alleles exist, but trans phase is not proven"
                    )
            else:
                for decision in decisions:
                    decision["decision_status"] = "SECOND_ALLELE_REQUIRED"

        rows: list[dict[str, Any]] = []
        for decision in self._decisions:
            decision = dict(decision)
            decision.pop("_transcripts", None)
            decision["decision_id"] = hashlib.sha256(canonical_json_bytes(decision)).hexdigest()
            rows.append(decision)
        rows.sort(key=lambda row: (row["variant_id"], row["gene_symbol"], row["policy_key"]))
        actionable = sorted({row["variant_id"] for row in rows if row["is_actionable"]})
        actionable_pairs = [
            {"variant_id": variant_id, "gene_symbol": gene_symbol}
            for variant_id, gene_symbol in sorted(
                {
                    (str(row["variant_id"]), str(row["gene_symbol"]))
                    for row in rows
                    if row["is_actionable"]
                }
            )
        ]
        status_counts: dict[str, int] = defaultdict(int)
        for row in rows:
            status_counts[str(row["decision_status"])] += 1
        return {
            "rows": rows,
            "actionable_variant_ids": actionable,
            "actionable_variant_gene_pairs": actionable_pairs,
            "status_counts": dict(sorted(status_counts.items())),
            "policy_version": self.version,
            "policy_snapshot": self.snapshot,
            "minimum_clinvar_review_stars": self.min_review_stars,
        }


def materialize_acmg_clinical_artifacts(
    evaluation: dict[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    """Map policy decisions into normative finding/evidence/interpretation rows."""

    findings: list[dict[str, Any]] = []
    policy_evidence: list[dict[str, Any]] = []
    interpretations: list[dict[str, Any]] = []
    for decision in evaluation["rows"]:
        decision_id = str(decision["decision_id"])
        finding_id = hashlib.sha256(f"finding\0{decision_id}".encode()).hexdigest()
        policy_evidence_id = hashlib.sha256(f"policy-evidence\0{decision_id}".encode()).hexdigest()
        interpretation_id = hashlib.sha256(f"interpretation\0{decision_id}".encode()).hexdigest()
        evidence_ids = sorted(
            [*{str(value) for value in decision["clinvar_evidence_ids"]}, policy_evidence_id]
        )
        findings.append(
            {
                "finding_id": finding_id,
                "variant_id": decision["variant_id"],
                "gene_symbol": decision["gene_symbol"],
                "condition": decision["disease"],
                "claim_type": "disease_risk",
                "classification": decision["clinvar_significance"],
                "clinical_grade": False,
                "evidence_ids": evidence_ids,
            }
        )
        policy_evidence.append(
            {
                "evidence_id": policy_evidence_id,
                "finding_id": finding_id,
                "variant_id": decision["variant_id"],
                "source": "ACMG_SF",
                "source_record_id": decision["policy_key"],
                "source_version": (
                    f"v{decision['policy_version']}/snapshot-{decision['policy_snapshot']}"
                ),
                "assertion": (
                    f"{decision['gene_symbol']} / {decision['disease']}: "
                    f"{decision['variants_to_report']}"
                ),
                "review_status": (
                    "published policy; application to this genotype is algorithmic and unreviewed"
                ),
                "retrieved_at": _snapshot_datetime(str(decision["policy_snapshot"])),
            }
        )
        genotype = decision["genotype"]
        separator = "|" if decision["phased"] else "/"
        if decision["decision_status"] == "REPORTABLE_CANDIDATE":
            genotype_matches_model: bool | None = True
        elif decision["decision_status"] in {
            "GENOTYPE_NOT_REPORTABLE",
            "SECOND_ALLELE_REQUIRED",
        }:
            genotype_matches_model = False
        else:
            genotype_matches_model = None
        interpretations.append(
            {
                "interpretation_id": interpretation_id,
                "finding_id": finding_id,
                "variant_id": decision["variant_id"],
                "gene_symbol": decision["gene_symbol"],
                "interpretation_status": "algorithmic",
                "claim_type": "disease_risk",
                "condition": decision["disease"],
                "inheritance_model": decision["inheritance"],
                "required_genotype": decision["variants_to_report"],
                "observed_genotype": separator.join(str(value) for value in genotype),
                "genotype_matches_model": genotype_matches_model,
                "penetrance": None,
                "effect_size": None,
                "analytical_confidence": decision["analytical_confidence"],
                "classification_review_status": (
                    f"ClinVar {decision['clinvar_review_stars']} star(s); "
                    f"conflicts={str(decision['clinvar_has_conflicts']).lower()}"
                ),
                "gene_disease_validity": None,
                "actionability_basis": (
                    f"ACMG SF v{decision['policy_version']} / {decision['decision_status']}"
                ),
                "confirmation_required": True,
                "reviewed_by": None,
                "reviewed_at": None,
                "evidence_snapshot": (
                    f"ClinVar {decision['clinvar_snapshot']}; "
                    f"ACMG/ClinGen {decision['policy_snapshot']}"
                ),
                "evidence_ids": evidence_ids,
            }
        )
    return {
        "clinical_findings": findings,
        "policy_evidence": policy_evidence,
        "interpretations": interpretations,
    }
