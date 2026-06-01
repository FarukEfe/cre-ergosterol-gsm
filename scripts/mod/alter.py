# Cobra package
from cobra import io
from cobra.core import Reaction, Metabolite, Model

# Other
import os, sys, argparse
import pandas as pd
from typing import Literal
from itertools import product

# Toolbox
from scripts.utils.parse_stoic import parse_stoichiometry

def alter(
    args, 
    ref: Model,
    ref_name: str,
    rxns_df: pd.DataFrame,
    genes_df: pd.DataFrame,
    sqe_source: Literal['Bb', 'Sc', 'Cf', None] = None, 
    sqs_source: Literal['Sc', 'Cl', 'Tv', None] = None,
    mva: bool = False
):

    # Safety
    if sqe_source not in ['Bb', 'Sc', 'Cf', None]:
        raise ValueError(f"Invalid SQE source '{sqe_source}'. Must be one of 'Bb', 'Sc', 'Cf', or None.")
    if sqs_source not in ['Sc', 'Cl', 'Tv', None]:
        raise ValueError(f"Invalid SQS source '{sqs_source}'. Must be one of 'Sc', 'Cl', 'Tv', or None.")

    # Add reactions and fill gene annotations (model already copied)
    print("Adding rxns ...\n")
    for _, row in rxns_df.iterrows():

        # Filter out if no in current construct
        if not (
            (sqe_source and (row['sbml_id'] == f'ALT_SQE{sqe_source}')) or
            (sqs_source and (
                row['sbml_id'] == f'ALT_SQS{sqs_source}' or
                row['sbml_id'] == f'ALT_PSPPS{sqs_source}'
                )
            ) or
            (mva and row['sbml_id'] in ['ALT_MVAS', 'ALT_MVAE', 'ALT_MVK', 'ALT_PMK', 'ALT_MVAD', 'ALT_IDLI'])
        ):
            continue

        try:
            # Ignore if reaction already in model (from sbml)
            rxn = ref.reactions.get_by_id(row['sbml_id'])
        except KeyError:
            # Add reaction fields to object
            rxn = Reaction(row['sbml_id'])
            rxn.name = row['display_name']
            rxn.lower_bound = float(row['lower_bound'])
            rxn.upper_bound = float(row['upper_bound'])
            rxn.subsystem = row['subsystem']
            rxn.annotation = {
                'kegg.reaction': row['kegg_reaction_id'],
                'ec-code': row['ec_number'],
            }
            # Add metabolites to reaction
            rxn.add_metabolites(parse_stoichiometry(eqn=row['stoichiometric_equation'], model=ref))
            # Add gpr to reaction if exists
            if pd.notna(row['gpr_rule']) and row['gpr_rule'].strip():
                rxn.gene_reaction_rule = row['gpr_rule']
            # Add reaction to model
            ref.add_reactions([rxn])

    print("Patching gene annotations ...\n")
    for _, row in genes_df.iterrows():

        try: 
            gene = ref.genes.get_by_id(row['sbml_id'])
            gene.name = row['display_name']
            gene.annotation = {
                'kegg.genes': row['annotation_kegg_genes'] if pd.notna(row['annotation_kegg_genes']) else '',
                'uniprot': row['uniprot_accession'] if pd.notna(row['uniprot_accession']) else '',
            }
        except KeyError: 
            continue # Gene not referenced by any reaction GPR. Skip.

    # Print out results
    new_name = f"{ref_name}{f'_SQS{sqs_source}' if sqs_source else ''}{f'_SQE{sqe_source}' if sqe_source else ''}{'_MVA' if mva else ''}"
    print(f"New model {new_name} has {len(ref.reactions)} reactions.")

    # Save altered model to repo
    save_path = os.path.join(args.outdir, ref_name)
    if not os.path.exists(save_path): os.makedirs(save_path)
    filepath = os.path.join(save_path, f"{new_name}.xml")
    io.write_sbml_model(ref, filepath)


if __name__ == "__main__":

    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='alter',
        description='Alter your model to specific constructs.'
    )
    parser.add_argument('--sbmlpath', required=True, help='Path to the input SBML model file.')
    parser.add_argument('--outdir', default="./data/altered/xmls", help='Directory to save the altered model SBML files.')
    parser.add_argument('--datadir', default="./data/altered/tables/stable", help='Directory containing the input CSV tables for reactions and compounds.')
    args = parser.parse_args()

    # Load sbml model
    ref, _ = io.validate_sbml_model(args.sbmlpath, validate=True)
    if not ref:
        print('No model recognized. Exiting...')
        sys.exit(1)

    # Extract model name
    ref_name = os.path.split(args.sbmlpath)[-1].split('.')[0]
    print(f'Old model {ref_name} has {len(ref.reactions)} reactions.')

    # Import alteration tables for reactions and compounds
    rxns_df = pd.read_csv(os.path.join(args.datadir, "reactions.csv"))
    mets_df = pd.read_csv(os.path.join(args.datadir, "metabolites.csv"))
    genes_df = pd.read_csv(os.path.join(args.datadir, "genes.csv"))

    # Add metabolites
    for _, row in mets_df.iterrows():
        try:
            met = ref.metabolites.get_by_id(row['sbml_id'])
        except KeyError:
            met = Metabolite(row['sbml_id'])
            met.name = row['display_name']
            met.formula = str(row['formula']) if pd.notna(row['formula']) else ''
            met.charge = int(row['charge']) if pd.notna(row['charge']) else 0
            met.compartment = row['compartment_id']
            met.annotation = {
                'kegg.compound': row['kegg_id'] if pd.notna(row['kegg_id']) else '',
                'chebi': row['chebi_id'] if pd.notna(row['chebi_id']) else '',
                'pubchem.compound': str(row['pubchem_cid']) if pd.notna(row['pubchem_cid']) else '',
                'inchikey': row['inchikey'] if pd.notna(row['inchikey']) else '',
            }
            ref.add_metabolites([met])

    # Alter
    sqe_sources = ['Bb', 'Sc', 'Cf', None]
    sqs_sources = ['Sc', 'Cl', 'Tv', None]
    sqe_sqs = list(product(sqe_sources, sqs_sources))
    
    for sqe_source, sqs_source in sqe_sqs:

        if (sqe_source, sqs_source) == (None, None):
            continue # Skip if no alteration, mva irrelevant w/o sqe or sqs

        # w/o mva
        alter(
            args, 
            ref=ref.copy(), ref_name=ref_name, 
            rxns_df=rxns_df, genes_df=genes_df, 
            sqe_source=sqe_source, sqs_source=sqs_source
        )

        # w/ mva
        alter(
            args, 
            ref=ref.copy(), ref_name=ref_name, 
            rxns_df=rxns_df, genes_df=genes_df, 
            sqe_source=sqe_source, sqs_source=sqs_source, mva=True
        )
