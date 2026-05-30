from cobra import io
from cobra.core import Metabolite, Reaction, Gene, Model
import os, argparse, pandas as pd

from scripts.utils.model import *
from scripts.utils.tools import split_coef, split_coef_reac
from scripts.utils.parse_stoic import parse_stoichiometry

def gapfill_model(args):

    old, err = io.validate_sbml_model(args.sbmlpath)
    if not old:
        print(f"Error loading model: {err}")
        exit(1)

    model = old.copy()
    # Get file name
    model_name = os.path.splitext(os.path.basename(args.sbmlpath))[0]

    # Import the manual reactions and compounds tables

    # sbml_id,display_name,phytozome_locus,kegg_gene_id,uniprot_accession,
    # ncbi_gene_id,kegg_ko,ec_number,gene_symbol,function,pathway_section,
    # annotation_kegg_genes,annotation_uniprot,annotation_ncbigene,notes
    genes_df = pd.read_csv(os.path.join(args.datadir, "genes.csv"))
    # sbml_id,display_name,stoichiometric_equation,lower_bound,upper_bound,
    # directionality,gpr_rule,ec_number,kegg_reaction_id,subsystem,
    # annotation_kegg_reaction,annotation_ec,pathway_section,
    # icre1355_equivalent_id,notes
    reactions_df = pd.read_csv(os.path.join(args.datadir, "reactions.csv"))
    # sbml_id,display_name,compartment_id,compartment_name,formula,charge,
    # kegg_id,chebi_id,pubchem_cid,cas_number,inchikey,molecular_weight,
    # pathway_section,annotation_kegg,annotation_chebi,annotation_pubchem,
    # annotation_inchikey,notes
    mets_df = pd.read_csv(os.path.join(args.datadir, "metabolites.csv"))

    print("Adding mets ...\n")
    for index, row in mets_df.iterrows():

        try:
            met = model.metabolites.get_by_id(row['sbml_id'])
        except KeyError:
            met = Metabolite(row['sbml_id'])
            met.name = row['display_name']
            met.formula = str(row['formula']) if pd.notna(row['formula']) else ''
            met.charge = int(row['charge']) if pd.notna(row['charge']) else 0
            met.compartment = row['compartment_id']
            met.annotation = {
                'kegg.compound': row['kegg_id'],
                'chebi': row['chebi_id'],
                'pubchem.compound': row['pubchem_cid'],
                'inchikey': row['inchikey'],
            }
            model.add_metabolites([met])

    print("Adding rxns ...\n")
    for index, row in reactions_df.iterrows():

        try:
            # Ignore if reaction already in model (from sbml)
            rxn = model.reactions.get_by_id(row['sbml_id'])
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
            rxn.add_metabolites(parse_stoichiometry(eqn=row['stoichiometric_equation'], model=model))
            # Add gpr to reaction if exists
            if pd.notna(row['gpr_rule']) and row['gpr_rule'].strip():
                rxn.gene_reaction_rule = row['gpr_rule']
            # Add reaction to model
            model.add_reactions([rxn])

    print("Patching gene annotations ...\n")
    for index, row in genes_df.iterrows():

        try: 
            gene = model.genes.get_by_id(row['sbml_id'])
            gene.name = row['display_name']
            gene.annotation = {
                'kegg.genes': row['annotation_kegg_genes'],
                'uniprot': row['uniprot_accession'],
            }
        except KeyError: 
            continue # Gene not referenced by any reaction GPR. Skip.

    # Debug
    print(f"\n\nFinal model has {len(model.metabolites)} metabolites and {len(model.reactions)} reactions.")
    print(f"Old model has {len(old.metabolites)} metabolites and {len(old.reactions)} reactions.")

    # Save model
    os.makedirs(args.outdir, exist_ok=True)
    save_file = os.path.join(args.outdir, f"MNL_{model_name}_GAPFILL.xml")
    io.write_sbml_model(model, save_file)


if __name__ == "__main__":
    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='_load_model',
        description='Load and validate your fba metabolic model from the .sbml format.'
    )
    parser.add_argument('--sbmlpath', default="./data/raw/iCre1355/iCre1355_auto.xml", help='Path to the input SBML model file.')
    parser.add_argument('--outdir', default="./data/fill/xmls", help='Directory to save the gapfilled model SBML file.')
    parser.add_argument('--datadir', default="./data/fill/tables/stable", help='Directory containing the input CSV tables for genes, reactions, and metabolites.')
    args = parser.parse_args()

    gapfill_model(args)