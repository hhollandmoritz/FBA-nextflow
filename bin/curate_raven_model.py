#!/usr/bin/env python3
"""Apply reproducible reaction changes to a COBRA model."""

import argparse
import csv

import cobra


def read_table(path):
    """Read a TSV file as dictionaries."""

    if path is None:
        return []

    with open(path, newline="") as handle:
        return list(csv.DictReader(
            handle,
            delimiter="\t",
        ))


def remove_reactions(model, rows, changes):
    """Remove specified reactions."""

    for row in rows:
        reaction_id = row["reaction_id"]

        if reaction_id not in model.reactions:
            raise ValueError(
                f"Reaction not found: {reaction_id}"
            )

        model.remove_reactions(
            [reaction_id],
            remove_orphans=False,
        )

        changes.append({
            "reaction_id": reaction_id,
            "operation": "remove",
            "details": "Removed from model",
        })


def modify_reactions(model, rows, changes):
    """Modify reaction bounds."""

    for row in rows:
        reaction_id = row["reaction_id"]

        if reaction_id not in model.reactions:
            raise ValueError(
                f"Reaction not found: {reaction_id}"
            )

        reaction = model.reactions.get_by_id(reaction_id)

        old_bounds = reaction.bounds
        lower, upper = old_bounds

        if row.get("lower_bound"):
            lower = float(row["lower_bound"])

        if row.get("upper_bound"):
            upper = float(row["upper_bound"])

        reaction.bounds = (lower, upper)

        changes.append({
            "reaction_id": reaction_id,
            "operation": "modify",
            "details": f"{old_bounds} -> {reaction.bounds}",
        })


def write_report(changes, output):
    """Save an audit trail of all modifications."""

    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "reaction_id",
                "operation",
                "details",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(changes)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--model", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--remove-reactions")
    parser.add_argument("--modify-reactions")
    parser.add_argument("--report", required=True)

    args = parser.parse_args()

    model = cobra.io.read_sbml_model(args.model)
    changes = []

    remove_reactions(
        model,
        read_table(args.remove_reactions),
        changes,
    )

    modify_reactions(
        model,
        read_table(args.modify_reactions),
        changes,
    )

    cobra.io.write_sbml_model(model, args.output)
    write_report(changes, args.report)

    print(f"Applied {len(changes)} changes")


if __name__ == "__main__":
    main()
