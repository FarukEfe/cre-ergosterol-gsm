from cobra.flux_analysis import flux_variability_analysis
from cobra.util import create_stoichiometric_matrix
from cobra.core import Model, Gene, Metabolite, Reaction

import numpy as np
from scipy.linalg import null_space
import networkx as nx

def gene_in_model(model: Model, gene: str):
    try:
        model.genes.get_by_id(gene)
        return True
    except:
        return False

def rxn_in_model(model: Model, rxn: str):
    try:
        model.reactions.get_by_id(rxn)
        return True
    except:
        return False

def met_in_model(model: Model, met: str):
    try:
        model.metabolites.get_by_id(met)
        return True
    except:
        return False

def add_single_gene_reaction_pair(
    model: Model, 
    gene_id: str, 
    reaction_id: str, 
    reaction_name: str,
    reaction_subsystem: str, 
    metabolites: list[tuple[int, str]],
    gene_name=None,
    reversible=False
):
    """Add a gene-reaction pair to the model."""

    # Avoid duplicates for gene and reaction
    # assert not model.genes.query(lambda k: k == gene_id, attribute='id')
    assert not model.reactions.query(lambda k: k == reaction_id, attribute='id')
    # Avoid metabolites not in model and subsystem exists
    assert all(list(map(lambda x: met_in_model(model, x[1]), metabolites)))
    assert reaction_subsystem in [grp.name for grp in model.groups]

    # Add gene and reaction to model
    rxn = Reaction(id=reaction_id)

    if gene_name is None:
        gene_name = gene_id
    
    # Set gene
    gene = None
    if gene_id in [g.id for g in model.genes]:
        gene = model.genes.get_by_id(gene_id)
    else:
        gene = Gene(gene_id, name=gene_name)
        model.genes.add(gene)

    # Set reaction
    model.add_reactions([rxn])

    rxn.name = reaction_name
    rxn.bounds = (-1000, 1000) if reversible else (0, 1000)

    # Find metabolite objects from id list
    add_mets = {}
    for coeff, met_id in metabolites:
        met = model.metabolites.get_by_id(met_id)
        add_mets[met] = coeff
    # Add metabolites and gene to reaction data
    rxn.add_metabolites(add_mets)
    rxn.gene_reaction_rule = gene_id
