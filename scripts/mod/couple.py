from cobra import io
from cobra.core import Metabolite, Reaction, Gene, Model
import os, argparse, sys

def couple_biomass(args):
    
    # Load Model
    ref, _ = io.validate_sbml_model(args.sbmlpath, validate=True)
    if not ref:
        raise ValueError('No model recognized. Exiting...')

    # Copy new instance
    new = ref.copy()

    # Couple biomass to ergosterol pathway
    biomass_auto = new.reactions.get_by_id('Biomass_Chlamy_auto')
    biomass_mixo = new.reactions.get_by_id('Biomass_Chlamy_mixo')
    biomass_hetero = new.reactions.get_by_id('Biomass_Chlamy_hetero')

    # Coefficients in mmol/gDW (from Voshall 2021, total 3.4 nmol/mg DW)
    # i.e. 0.0034 mmol/gDW total ergosterol. 
    # Assume that pool is demanded at 100% in all biomass modes.
    sterol_coefficients = {
        'ergosterol_c':          -0.002142,   # 63%
        '7dhporiferasterol_c':   -0.000918,   # 27%
        'ergost7enol_c':         -0.000340,   # 10%
    }

    # Couple ergosterol
    biomass_auto.add_metabolites(sterol_coefficients)
    biomass_mixo.add_metabolites(sterol_coefficients)
    biomass_hetero.add_metabolites(sterol_coefficients)

    # Save new model
    os.makedirs(args.outdir, exist_ok=True)
    save_path = os.path.join(args.outdir, f"{os.path.basename(args.sbmlpath)}")
    io.write_sbml_model(new, save_path)
    print(f"New model with coupled biomass saved to {save_path}.")


if __name__ == "__main__":

    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='couple',
        description='Couple the biomass reaction to the ergosterol pathway.'
    )
    parser.add_argument('--sbmlpath', default="./data/fill/xmls/MNL_iCre1355_auto_GAPFILL.xml", help='Path to the input SBML model file.')
    parser.add_argument('--outdir', default="./data/coupled/xmls", help='Directory to save the coupled model SBML file.')
    args = parser.parse_args()

    try:
        couple_biomass(args)
    except ValueError as error:
        print(error)
        sys.exit(1)
