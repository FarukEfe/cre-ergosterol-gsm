from cobra import io
from cobra import Metabolite, Reaction, Gene, Model
import os, argparse, sys

def sink(args):

    # Create directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    
    # Load model
    ref, _ = io.validate_sbml_model(args.sbmlpath, validate=True)
    if not ref:
        raise ValueError('No model recognized. Exiting...')
    
    new = ref.copy()
    for met_id in args.metabolites:
        try:
            met = new.metabolites.get_by_id(met_id)
            sink_rxn = Reaction(f'sink_{met_id}')
            sink_rxn.name = f'Sink for {met_id}'
            sink_rxn.lower_bound = 0
            sink_rxn.upper_bound = 1000
            sink_rxn.add_metabolites({met: -1})
            new.add_reactions([sink_rxn])
        except KeyError:
            raise ValueError(f"Metabolite {met_id} not found in model. Exiting...")
            sys.exit(1)

    # Save new model
    os.makedirs(args.outdir, exist_ok=True)
    save_path = os.path.join(args.outdir, f"{os.path.basename(args.sbmlpath)}")
    io.write_sbml_model(new, save_path)
    print(f"New model with sink reactions saved to {save_path}.")

if __name__ == "__main__":

    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='sink',
        description='Create a sink/demand reaction for given metabolite(s) in the model.'
    )
    parser.add_argument('--sbmlpath', default="./data/coupled/xmls/MNL_iCre1355_auto_GAPFILL.xml", help='Path to the input SBML model file.')
    parser.add_argument('--metabolites', nargs='+', default=['ergosterol_c', '7dhporiferasterol_c', 'ergost7enol_c'], help='List of metabolite IDs to create sink reactions for.')
    parser.add_argument('--outdir', default="./data/sink/xmls", help='Directory to save the new model SBML file.')
    args = parser.parse_args()

    try:
        sink(args)
    except ValueError as error:
        print(error)
        sys.exit(1)