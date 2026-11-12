#!/usr/bin/env python3
"""Add, modify, and remove reactions from a COBRA-compatible SBML model.

Steps:
  1. Add metabolites                  (--metabolites)
  2. Import reference-model reactions (--import-reactions, --reference-model)
  3. Add custom reactions             (--add-reactions)
  4. Modify bounds and/or GPRs        (--modify-reactions)
  5. Remove reactions                 (--remove-reactions)
  6. Set an objective (optional)      (--objective)

"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

import cobra

# Columns and files needed for creating the report of changes
AUDIT_COLUMNS = ["operation", "reaction_id", "before", "after", "source", "evidence"]
 
def record(changes, operation, reaction_id, before="", after="", source="", evidence=""):
    """Utility function that records a change to the model"""
    changes.append({
        "operation": operation, "reaction_id": reaction_id,
        "before": str(before), "after": str(after),
        "source": source or "", "evidence": evidence or "",
    })

def add_metabolites(model, rows, changes):
    """Add metabolites defined by ID, compartment and optional metadata."""
    for row in rows:
        mid = (row.get("metabolite_id") or "").strip()
        compartment = (row.get("compartment") or "").strip()
        if not mid or not compartment:
            raise ValueError("Metabolite rows require metabolite_id and compartment")
        if mid in model.metabolites:
            raise ValueError(f"Metabolite already exists: {mid}")
        kwargs = {"name": (row.get("name") or "").strip() or mid,
                  "compartment": compartment}
        formula = (row.get("formula") or "").strip()
        if formula:
            kwargs["formula"] = formula
        charge_text = (row.get("charge") or "").strip()
        if charge_text:
            try:
                kwargs["charge"] = int(charge_text)
            except ValueError as exc:
                raise ValueError(f"{mid}: charge must be an integer") from exc
        met = cobra.Metabolite(mid, **kwargs)
        model.add_metabolites([met])
        record(changes, "add_metabolite", mid, after=compartment,
               source=row.get("source"), evidence=row.get("evidence"))


def import_reactions(model, reference, rows, changes):
    """Copy reactions from a reference SBML using a list provided by the user;
    Reaction metabolites are matched by ID and compartment.
    """
    if rows and reference is None:
        raise ValueError("--import-reactions requires --reference-model")
    for row in rows:
        rid = (row.get("reaction_id") or "").strip()
        if not rid:
            raise ValueError("Import reaction row has an empty reaction_id")
        if rid in model.reactions:
            raise ValueError(f"Did not add reaction. Reaction already exists: {rid}")
        if rid not in reference.reactions:
            raise ValueError(f"Reaction not found in reference SBML: {rid}")
        original = reference.reactions.get_by_id(rid)
        stoich = {}
        for met, coeff in original.metabolites.items():
            if met.id in model.metabolites:
                target = model.metabolites.get_by_id(met.id)
                if target.compartment != met.compartment:
                    raise ValueError(
                        f"{rid}: metabolite {met.id} compartment differs "
                        f"({target.compartment} versus {met.compartment})"
                    )
            else:
                target = met.copy()
                model.add_metabolites([target])
            stoich[target] = coeff
        rxn = cobra.Reaction(rid, name=original.name,
                             lower_bound=original.lower_bound,
                             upper_bound=original.upper_bound)
        rxn.add_metabolites(stoich)
        rxn.gene_reaction_rule = original.gene_reaction_rule
        rxn.annotation = dict(original.annotation)
        rxn.notes = dict(original.notes)
        model.add_reactions([rxn])
        record(changes, "import_reaction", rid, after=rxn.reaction,
               source=row.get("source") or "reference_model",
               evidence=row.get("evidence"))


def add_reactions(model, rows, changes):
    """Add custom reactions using a JSON {metabolite_id: coefficient} map."""
    for row in rows:
        rid = (row.get("reaction_id") or "").strip()
        if not rid:
            raise ValueError("Custom reaction row has an empty reaction_id")
        if rid in model.reactions:
            raise ValueError(f"Did not add reaction. Reaction already exists: {rid}")
        try:
            raw = json.loads(row["stoichiometry"])
        except (ValueError, TypeError, KeyError) as exc:
            raise ValueError(f"{rid}: stoichiometry must be a JSON object") from exc
        if not isinstance(raw, dict) or not raw:
            raise ValueError(f"{rid}: stoichiometry must be a nonempty JSON object")
        unknown = sorted(mid for mid in raw if mid not in model.metabolites)
        if unknown:
            raise ValueError(f"{rid}: undefined metabolites: {unknown}")
        coefficients = {
            model.metabolites.get_by_id(mid): float(value)
            for mid, value in raw.items()
        }
        if not any(x < 0 for x in coefficients.values()) or not any(x > 0 for x in coefficients.values()):
            raise ValueError(f"{rid}: reaction must have reactants and products")
        lower = float(row.get("lower_bound") or 0)
        upper = float(row.get("upper_bound") or 1000)
        if lower > upper:
            raise ValueError(f"{rid}: lower_bound exceeds upper_bound")
        rxn = cobra.Reaction(rid, name=(row.get("name") or "").strip() or rid,
                             lower_bound=lower, upper_bound=upper)
        rxn.add_metabolites(coefficients)
        rxn.gene_reaction_rule = (row.get("gpr") or "").strip()
        model.add_reactions([rxn])
        record(changes, "add_reaction", rid, after=rxn.reaction,
               source=row.get("source"), evidence=row.get("evidence"))


def modify_reactions(model, rows, changes):
    """Change reaction bounds and optionally its GPR; blank values mean unchanged."""
    for row in rows:
        rid = (row.get("reaction_id") or "").strip()
        if rid not in model.reactions:
            raise ValueError(f"Can't modify. Reaction not found: {rid}")
        rxn = model.reactions.get_by_id(rid)
        before = f"bounds={rxn.bounds}; gpr={rxn.gene_reaction_rule}"
        lower, upper = rxn.bounds
        if (row.get("lower_bound") or "").strip():
            lower = float(row["lower_bound"])
        if (row.get("upper_bound") or "").strip():
            upper = float(row["upper_bound"])
        if lower > upper:
            raise ValueError(f"{rid}: lower_bound exceeds upper_bound")
        rxn.bounds = (lower, upper)
        if (row.get("gpr") or "").strip():
            rxn.gene_reaction_rule = row["gpr"].strip()
        after = f"bounds={rxn.bounds}; gpr={rxn.gene_reaction_rule}"
        record(changes, "modify_reaction", rid, before=before, after=after,
               source=row.get("source"), evidence=row.get("evidence"))


def remove_reactions(model, rows, changes):
    """Remove specified reactions. (Keeps orphan metabolites)."""
    for row in rows:
        rid = (row.get("reaction_id") or "").strip()
        if rid not in model.reactions:
            raise ValueError(f"Can't remove. Reaction not found: {rid}")
        rxn = model.reactions.get_by_id(rid)
        before = rxn.reaction
        model.remove_reactions([rxn], remove_orphans=False)
        record(changes, "remove_reaction", rid, before=before,
               source=row.get("source"), evidence=row.get("evidence"))


def write_report(changes, output):
    """Save one audit row for each modification."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(changes, columns=AUDIT_COLUMNS).to_csv(output, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True, help="Draft SBML model")
    parser.add_argument("--output", type=Path, required=True, help="Adapted SBML model")
    parser.add_argument("--report", type=Path, required=True, help="Curation audit TSV")
    parser.add_argument("--metabolites", type=Path, help="TSV of custom metabolites")
    parser.add_argument("--import-reactions", type=Path, help="TSV of reaction IDs from reference to add or modify")
    parser.add_argument("--reference-model", type=Path, help="SBML source for imported reactions")
    parser.add_argument("--add-reactions", type=Path, help="TSV of custom reactions")
    parser.add_argument("--modify-reactions", type=Path, help="TSV of bounds/GPR edits")
    parser.add_argument("--remove-reactions", type=Path, help="TSV of reaction IDs to delete")
    parser.add_argument("--objective", help="Optional reaction ID to set as the objective")
    args = parser.parse_args()

    model = cobra.io.read_sbml_model(str(args.model))
    reference = (cobra.io.read_sbml_model(str(args.reference_model))
                 if args.reference_model is not None else None)
    
    # initialize a list to track changes made to the model
    changes = []

    # Load each TSV and extract the rows for handing off to the modificaiton funcitons.
    required_columns = {
        "metabolites": ("metabolite_id", "compartment"),
        "add_reactions": ("reaction_id", "stoichiometry"),
    }
    tables = {}
    for key in ("metabolites", "import_reactions", "add_reactions",
                "modify_reactions", "remove_reactions"):
        path = getattr(args, key)
        # If no file is provided for a particular modifcation type, we'll just
        # use a blank list (so that modifications are optional).
        if path is None:
            tables[key] = []
            continue

        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                         encoding="utf-8-sig")
        required = required_columns.get(key, ("reaction_id",))
        missing = [col for col in required if col not in df.columns]
        if missing:
            raise ValueError(f"{path}: missing required TSV columns {missing}")

        # Ignore blank lines.
        df = df.loc[df.apply(lambda col: col.str.strip()).ne("").any(axis=1)]
        tables[key] = df.to_dict("records")
    # Call each type of modification function on the model
    add_metabolites(model, tables["metabolites"], changes)
    import_reactions(model, reference, tables["import_reactions"], changes)
    add_reactions(model, tables["add_reactions"], changes)
    modify_reactions(model, tables["modify_reactions"], changes)
    remove_reactions(model, tables["remove_reactions"], changes)

    # Modify the objective reaction
    if args.objective is not None:
        if args.objective not in model.reactions:
            raise ValueError(f"User-specified objective reaction not found in model: {args.objective}")
        before = str(model.objective.expression)
        model.objective = args.objective
        record(changes, "set_objective", args.objective, before=before,
               after=args.objective)

    # Write outputs: Create output directory, write modified model and changes report.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(model, str(args.output))
    write_report(changes, args.report)
    print(f"Applied {len(changes)} edits to {model.id}; wrote {args.output}")


if __name__ == "__main__":
    main()
