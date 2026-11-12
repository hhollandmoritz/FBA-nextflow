#!/usr/bin/env python3
"""Run RAVEN fast-LP gap-filling on an adapted SBML model.
Uses a template to fill gaps
"""

import argparse
from pathlib import Path

import cobra
import pandas as pd
from raven_toolbox.gapfilling import fill_gaps_fast_lp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, required=True, help='Adapted SBML model')
    parser.add_argument('--template', type=Path, required=True, help='SBML reaction template')
    parser.add_argument('--output', type=Path, required=True, help='Gap-filled SBML model')
    parser.add_argument('--report', type=Path, required=True, help='Added-reaction report (TSV)')
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()

    model = cobra.io.read_sbml_model(str(args.model))
    template = cobra.io.read_sbml_model(str(args.template))

    result = fill_gaps_fast_lp(model, template, n_proc=args.threads)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    cobra.io.write_sbml_model(result.model, str(args.output))

    args.report.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({'reaction_id': sorted(result.added_reactions)}).to_csv(
        args.report, sep='\t', index=False
    )
    print(f'Added {len(result.added_reactions)} reactions; wrote {args.output}')


if __name__ == '__main__':
    main()
