"""
This module is to parse the stoichiometry equation of a reaction to separate metabolites with their coefficients.
"""
import re, os, sys, argparse
from cobra import io, Model, Metabolite, Reaction


def parse_stoichiometry(eqn: str, model: Model, debug: bool = False) -> dict[Metabolite, float]:

    metabolites = {}

    reactants, products = eqn.split('-->')
    reactants = reactants.strip().split('+')
    products = products.strip().split('+')

    for token in reactants:
        token = token.strip()
        rmatch = re.match(r'^\((\d+\.?\d*)\)\s+(.+)$', token)
        if rmatch:
            coef = float(rmatch.group(1))
            met = rmatch.group(2)
        else:
            coef = 1.0
            met = token.strip()

        try:
            met = model.metabolites.get_by_id(met)
        except KeyError:
            raise ValueError(f"Metabolite '{met}' not found in model.")
        
        metabolites[met] = -1 * coef

        if debug: print(f"Reactant: {met.id}, Coefficient: {-1*coef}")

    for token in products:
        token = token.strip()
        pmatch = re.match(r'^\((\d+\.?\d*)\)\s+(.+)$', token)
        if pmatch:
            coef = float(pmatch.group(1))
            met = pmatch.group(2)
        else:
            coef = 1.0
            met = token.strip()

        try:
            met = model.metabolites.get_by_id(met)
        except KeyError:
            raise ValueError(f"Metabolite '{met}' not found in model.")

        metabolites[met] = coef

        if debug: print(f"Product: {met.id}, Coefficient: {coef}")

    return metabolites

# i.e. "(2) A + (3) B --> (4) C + (5) D" --> {"A": 2, "B": 3, "C": 4, "D": 5}
if __name__ == "__main__":
    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='parse_stoic',
        description='Parse the stoichiometry equation of a reaction to separate metabolites with their coefficients.'
    )
    parser.add_argument('--eqn', required=True, help='The stoichiometry equation to parse.')
    parser.add_argument('--model', required=True, help='The path to the SBML model file.')
    args = parser.parse_args()

    model = io.read_sbml_model(args.model)
    print(parse_stoichiometry(args.eqn, model))