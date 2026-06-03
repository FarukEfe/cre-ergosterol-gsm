# FVA
    # Run FVA
    # List blocked reactions
    # Export flux ranges
from cobra import io
import argparse, os
import warnings, sys

from cobra.core import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis

def run_flux_variability_analysis(
        model: Model,
        objective: str,
        loopless: bool = True,
        fraction_of_optimum: float = 0.9,
        reactions: list[str] = None
    ):
    """
    Perform flux variability analysis on the given model.

    Args:
        model (Model): The metabolic model to analyze.
        objective (str): The objective reaction for the analysis.
        reactions (list[str], optional): A list of reaction IDs to include in the analysis.

    Returns:
        pd.DataFrame: A DataFrame containing the flux ranges for each reaction.
    """

    model.objective = objective

    # cobra.flux_analysis.flux_variability_analysis
    solution = flux_variability_analysis(
        model,
        loopless=loopless,
        fraction_of_optimum=fraction_of_optimum,
        reaction_list=reactions,
        processes=None
    )
    return solution

if __name__ == "__main__":

    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='_load_model',
        description='Load and validate your fba metabolic model from the .sbml format.'
    )
    parser.add_argument('--sbmlpath', required=True, help='Path to the input SBML model file.')
    parser.add_argument('--objective', required=True, help='Objective reaction(s) for FVA, comma-separated.')
    parser.add_argument('-o', '--outdir', default="./res/fva", help='Directory to save the FVA results.')
    parser.add_argument('-l', '--loopless', action='store_true', help='Use loopless FVA.')
    parser.add_argument('-r', '--reactions', help='Reactions to include in the analysis, comma-separated.')
    args = parser.parse_args()

    print("Model import from SBML file... {}".format(args.sbmlpath))

    # Model Import
    model, error = io.validate_sbml_model(args.sbmlpath)

    if not model:
        print('No model recognized. Exiting...')
        sys.exit(1)

    # r: Reaction = model.reactions.get_by_id("SS")
    # r.objective_coefficient = 1.0
    flux_ranges = None
    reactions = args.reactions.split(',') if args.reactions else None

    try:
        flux_ranges = run_flux_variability_analysis(
            model,
            loopless=args.loopless,
            objective=args.objective,
            reactions=reactions
        )
    except Exception as e:
        print(f"Error during flux variability analysis: {e}")
        sys.exit(1)

    file_name: str = os.path.split(args.sbmlpath)[-1].split('.')[0]
    output_dest: str = os.path.join(args.outdir, file_name)
    if not os.path.exists(output_dest): os.makedirs(output_dest)

    export_path = os.path.join(output_dest, f"{args.objective}_fva.csv")
    flux_ranges.to_csv(export_path)
    exit(0)
