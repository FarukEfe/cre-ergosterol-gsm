from cobra import Reaction, Model, Metabolite, io
import numpy as np
from builtins import map

def split_coef(inp: str) -> tuple[int, str]:
    """Split coefficient from a product in equation (built for the repo .csv format)"""
    try:
        coef, cpd = inp[:1], inp[1:]
        coef = 1 if coef == "" else int(coef)
        return coef, cpd
    except:
        return 1, inp

def split_coef_reac(inp: str) -> tuple[int, str]:
    """Split coefficient from a reaction in equation (built for the repo .csv format)"""
    res = split_coef(inp)
    return -1 * res[0], res[1]

def find_rxns_with_metabolites(model: Model, metabolites: list[str]):

    rxns_with_mets = []
    for rxn in model.reactions:
        for met, coef in list(rxn.metabolites.items()):
            print(f'Rxn: {rxn.id} --- Met: {met.id}')
            if met.id in metabolites:
                rxns_with_mets.append(rxn.id)
                break
    return rxns_with_mets

def get_rxn_metabolites(reaction_id: str, model: Model) -> tuple[list[Metabolite], list[int]]:

    rxn = model.reactions.get_by_id(reaction_id)

    mets, coefs = [], []
    for met, coef in rxn.metabolites.items():
        if coef != 0:
            mets.append(met)
            coefs.append(coef)
    return mets, coefs