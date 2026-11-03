#!/usr/bin/env python3
"""Construct a KEGG-based RAVEN draft from protein sequences or precomputed KOs.

Requires raven-toolbox (RAVEN 3 Python), cobra, and pandas.
The KO route requires the *matching* KEGG reference SBML + KO-reaction table
(and optionally reaction quality flags) from the same RAVEN data release.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from importlib.metadata import version

import cobra
import pandas as pd

from raven_toolbox.reconstruction.kegg import (
    get_kegg_model_from_sequences_with_artefacts,
)

from raven_toolbox.reconstruction.kegg.assemble import (
    assemble_model_from_ko_genes,
)

KO_PATTERN = re.compile(r"(?<![A-Za-z0-9])K\d{5}(?!\d)")


def read_ko_assignments(path: Path, gene_column: str, ko_column: str):
    """Read a DRAM or normalized TSV, allowing multiple KO IDs in each cell."""
    assignments = defaultdict(set)
    total_rows = 0
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or gene_column not in reader.fieldnames or ko_column not in reader.fieldnames:
            raise ValueError(
                f"Expected columns {gene_column!r} and {ko_column!r}; "
                f"found {reader.fieldnames!r} in {path}"
            )
        for row_no, row in enumerate(reader, start=2):
            total_rows += 1
            gene = (row.get(gene_column) or "").strip()
            ko_text = (row.get(ko_column) or "").strip()
            if not gene or not ko_text or ko_text.lower() in {"na", "nan", "none", "null"}:
                continue
            kos = KO_PATTERN.findall(ko_text)
            if not kos:
                raise ValueError(
                    f"Row {row_no} contains a nonempty KO field with no recognizable "
                    f"Kxxxxx identifier: {ko_text!r}; check --ko-column"
                )
            for ko in kos:
                assignments[ko].add(gene)
    return {k: sorted(v) for k, v in sorted(assignments.items())}, total_rows


def read_table(path: Path):
    return pd.read_csv(path, sep="\t", dtype=str, compression="infer")


def write_evidence(model, output: Path):
    fields = ["reaction_id", "reaction_name", "equation", "gpr", "gene_ids",
              "lower_bound", "upper_bound", "annotation_json"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(fields)
        for rxn in sorted(model.reactions, key=lambda r: r.id):
            writer.writerow([
                rxn.id, rxn.name, rxn.reaction, rxn.gene_reaction_rule,
                ";".join(sorted(g.id for g in rxn.genes)),
                rxn.lower_bound, rxn.upper_bound,
                json.dumps(rxn.annotation, sort_keys=True, default=str),
            ])


def write_ko_evidence(ko_to_genes, output):
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["ko_id", "gene_id"])
        for ko, genes in sorted(ko_to_genes.items()):
            for gene in genes:
                writer.writerow([ko, gene])


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path,
                        help="Protein FASTA (required for sequence mode and auto fallback)")
    parser.add_argument("--ko-annotations", type=Path,
                        help="Optional DRAM or normalized gene-to-KO TSV")
    parser.add_argument("--gene-column", default="gene_id")
    parser.add_argument("--ko-column", default="ko_id")
    parser.add_argument("--method", choices=("auto", "ko", "sequence"), default="auto")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--ko-evidence", type=Path,
                        help="Optional normalized KO-to-gene assignment table")
    parser.add_argument("--model-id")
    parser.add_argument("--artefact-dir", type=Path,
                        help="RAVEN cached KEGG artefacts; recommended for offline HPC")
    parser.add_argument("--library", type=Path,
                        help="RAVEN HMM library for sequence reconstruction")
    parser.add_argument("--reference-model", type=Path,
                        help="KEGG reference SBML from same artefact release; required for KO route")
    parser.add_argument("--ko-reaction", type=Path,
                        help="KO-to-reaction TSV from same artefact release; required for KO route")
    parser.add_argument("--rxn-flags", type=Path,
                        help="Optional quality flag TSV from same artefact release")
    parser.add_argument("--domain", choices=("prokaryotes", "eukaryotes"),
                        default="prokaryotes")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--cutoff", type=float, default=1e-30)
    parser.add_argument("--min-score-ratio-ko", type=float, default=0.3)
    parser.add_argument("--min-score-ratio-g", type=float, default=0.9)
    parser.add_argument("--keep-spontaneous", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--keep-undefined-stoich", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--keep-incomplete", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--keep-general", action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    return args


def main():
    args = parse_args()
    if args.fasta is not None and not args.fasta.is_file():
        raise FileNotFoundError(args.fasta)
    ko_to_genes, n_rows = {}, 0

    if args.method != "sequence" and args.ko_annotations is not None:
        if args.ko_annotations.is_file():
            ko_to_genes, n_rows = read_ko_assignments(
                args.ko_annotations, args.gene_column, args.ko_column
            )
        elif args.method == "ko":
            raise FileNotFoundError(args.ko_annotations)
        # In auto mode a missing file is equivalent to no annotations.

    if args.method == "ko" and not ko_to_genes:
        raise ValueError("KO mode requires a TSV containing at least one valid KO assignment")

    selected = "ko" if args.method != "sequence" and ko_to_genes else "sequence"
    if selected == "sequence" and args.fasta is None:
        raise ValueError("Sequence reconstruction requires --fasta; supply a FASTA for auto fallback")
    model_id = args.model_id or (args.fasta.stem if args.fasta else args.ko_annotations.stem)
    common = dict(
        keep_spontaneous=args.keep_spontaneous,
        keep_undefined_stoich=args.keep_undefined_stoich,
        keep_incomplete=args.keep_incomplete,
        keep_general=args.keep_general,
    )

    if selected == "ko":
        if not args.reference_model or not args.ko_reaction:
            raise ValueError(
                "KO reconstruction requires --reference-model and --ko-reaction "
                "from the same RAVEN KEGG artefact release"
            )
        reference = cobra.io.read_sbml_model(str(args.reference_model))
        ko_reaction = read_table(args.ko_reaction)
        flags = read_table(args.rxn_flags) if args.rxn_flags else None
        model, _gpr_map = assemble_model_from_ko_genes(
            reference_model=reference,
            ko_reaction=ko_reaction,
            ko_to_genes=ko_to_genes,
            rxn_flags=flags,
            model_id=model_id,
            prune_orthology=True,
            **common,
        )
    else:
        kwargs = dict(
            domain=args.domain,
            threads=args.threads,
            cutoff=args.cutoff,
            min_score_ratio_ko=args.min_score_ratio_ko,
            min_score_ratio_g=args.min_score_ratio_g,
            model_id=model_id,
            **common,
        )
        model = get_kegg_model_from_sequences_with_artefacts(
            fasta=args.fasta,
            artefact_dir=args.artefact_dir,
            library=args.library,
            **kwargs,
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.evidence.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.parent.mkdir(parents=True, exist_ok=True)

    cobra.io.write_sbml_model(model, str(args.output))
    # Fail early if the model cannot round-trip through SBML.
    reloaded = cobra.io.read_sbml_model(str(args.output))
    if (len(reloaded.reactions), len(reloaded.metabolites), len(reloaded.genes)) != (
        len(model.reactions), len(model.metabolites), len(model.genes)
    ):
        raise RuntimeError("SBML round-trip changed model dimensions")
    write_evidence(model, args.evidence)
    if args.ko_evidence:
        args.ko_evidence.parent.mkdir(parents=True, exist_ok=True)
        write_ko_evidence(ko_to_genes if selected == "ko" else {}, args.ko_evidence)

    metadata = {
        "requested_method": args.method,
        "selected_method": selected,
        "fasta": str(args.fasta) if args.fasta else None,
        "ko_annotations": str(args.ko_annotations) if args.ko_annotations else None,
        "ko_input_rows": n_rows,
        "input_unique_kos": len(ko_to_genes) if selected == "ko" else None,
        "input_unique_ko_genes": len(set(g for genes in ko_to_genes.values() for g in genes)) if selected == "ko" else None,
        "model_id": model.id,
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "genes": len(model.genes),
        "raven_toolbox_version": version("raven-toolbox"),
        "cobra_version": version("cobra"),
        "artefact_dir": str(args.artefact_dir) if args.artefact_dir else None,
        "reference_model": str(args.reference_model) if selected == "ko" else None,
        "ko_reaction": str(args.ko_reaction) if selected == "ko" else None,
    }
    args.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
