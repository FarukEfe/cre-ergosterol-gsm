from cobra import io, Model, Solution
from cobra.flux_analysis import pfba
import escher
import argparse, os, sys

from cobra.util.solver import linear_reaction_coefficients

def run_flux_balance_analysis(
    model: Model,
    objective: str,
    is_pfba: bool = False,
    minimize: bool = False,
    fraction_of_optimum: float = 1.0,

):
    """
    Run Flux-Balance Analysis on the model for the provided objectives.
    Makes the assumption that all stated objectives in `objectives: list[str]` are equally important.

    Args:
        model (Model): The metabolic model to analyze.
        loopless (bool): Whether to use loopless FBA.
        pfba (bool): Whether to use parsimonious FBA.
        objective (str): The objective reaction for the analysis.
        reactions (list[str], optional): The list of reactions to include in the analysis.
        pfba_factor (float, optional): The factor to use for pfba. Defaults to 1.1.
        fraction_of_optimum (float, optional): The fraction of the optimum to use. Defaults to 1.0.

    Returns:
        Solution: The optimized solution for the provided objectives. For more see `cobra.Solution`
    """

    if is_pfba:
        model.objective = objective
        solution = pfba(
            model,
            objective={model.reactions.get_by_id(objective): 1.0},
            fraction_of_optimum=fraction_of_optimum,
        )
    else:
        model.objective = objective
        solution = model.optimize(raise_error=True)

    return solution


if __name__ == "__main__":
    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='_load_model',
        description='Load and validate your fba metabolic model from the .sbml format.'
    )
    parser.add_argument('--sbmlpath', required=True, help='Path to the input SBML model file.')
    parser.add_argument('--outdir', default="./res/fba", help='Directory to save the FBA results.')
    parser.add_argument('-p', '--pfba', action='store_true', help='Use parsimonious FBA.')
    parser.add_argument('--objective', required=True, help='Objective reaction for FBA.')
    args = parser.parse_args()

    # # Print arguments
    # print(f"Loading model from: {args.sbmlpath}")
    # print(f"Destination directory: {args.dest}")

    # Model Import
    model, error = io.validate_sbml_model(args.sbmlpath)

    if not model:
        print(f'Error loading model: {error}')
        sys.exit(1)
    
    solution = None
    try:
        # Run flux-balance analysis
        solution = run_flux_balance_analysis(model, objective=args.objective, is_pfba=args.pfba)
    except Exception as e:
        print(f"Error during FBA: {e}")
        sys.exit(1)

    # File name & destination
    file_name: str = os.path.split(args.sbmlpath)[-1].split('.')[0]
    dest_final: str = os.path.join(args.outdir, file_name)

    # Save to destination
    if not os.path.exists(dest_final): os.makedirs(dest_final)
    export_path = os.path.join(dest_final, f'{args.objective}.csv')
    solution.fluxes.to_csv(export_path)
    exit(0)