import argparse, json, os, sys
import numpy as np
import pandas as pd
from typing import TypedDict, Literal

from cobra import io, Model
from scripts.opt._fba import run_flux_balance_analysis
from scripts.opt._fva import run_flux_variability_analysis

class GrowthExperimentItem(TypedDict):
    model: Model
    mu_obs: float
    objective: str
    nh4_ub: float
    pi_ub: float
    acet_ub: float
    co2_ub: float
    light_ub: float | None

def predict_growth(
    model: Model, mu_obs: float, objective: str, 
    nh4_ub: float, pi_ub: float, acet_ub: float | Literal['NA'],
    co2_ub: float | Literal['ND'], light_ub: float | None = None
) -> tuple[float, float]:
    with model:

        if acet_ub == 'NA':
            model.reactions.get_by_id('EX_ac_e').bounds = (0.0, 0.0)
        else:
            model.reactions.get_by_id('EX_ac_e').lower_bound = acet_ub

        if co2_ub != 'ND':
            model.reactions.get_by_id('CO2t').upper_bound = co2_ub

        model.reactions.get_by_id('EX_nh4_e').lower_bound = nh4_ub
        model.reactions.get_by_id('EX_pi_e').lower_bound = pi_ub

        if light_ub is not None:
            model.reactions.get_by_id('PRISM_design_growth').upper_bound = light_ub

        model.objective = objective
        solution = model.optimize(raise_error=True)
        return mu_obs, solution.objective_value

def growth_experiment(
        experiments: dict[str, GrowthExperimentItem], 
        outdir: str
) -> tuple[pd.DataFrame, float]:
    items, mu_obs_l, mu_pred_l = [], [], []
    for name, item in experiments.items():
        print(f"Running experiment {name}...")
        mu_obs, mu_pred = predict_growth(**item)
        print(mu_obs, mu_pred)
        mu_obs_l.append(mu_obs)
        mu_pred_l.append(mu_pred)
        items.append({
            "name": name,
            "growth_mode": name.split('_')[0],
            "mu_obs": mu_obs,
            "mu_pred": mu_pred,
            # Exclude model from item dict to ensure .csv conversion
            **{k: v for k, v in item.items() if k != 'model'},
        })
    exit(1)
    # Conver to np.array
    mu_obs_l = np.array(mu_obs_l)
    mu_pred_l = np.array(mu_pred_l)
    # Compute R2
    mu_obs_mean = np.sum(mu_obs_l) / len(mu_obs_l)
    R2 = 1 - np.sum((mu_obs_l - mu_pred_l)**2)/np.sum((mu_obs_l - mu_obs_mean)**2)
    # Build dataframe and add to it
    df = pd.DataFrame(items)
    df["R2"] = R2
    # Save to .csv
    df.to_csv(os.path.join(outdir, "growth_experiment.csv"), index=False)
    # Return results 
    return df, R2
        

if __name__ == "__main__":

    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='_growth_experiment',
        description='Perform growth experiment on the model for the provided objectives.'
    )
    # Models of the Trophic Modes
    parser.add_argument('--model_type', default='sink', choices=['sink', 'coupled', 'fill', 'raw'], help='Model type to use for the growth experiment. Determines the SBML files to load.')
    # parser.add_argument('--sbml_auto', default='./data/sink/xmls/MNL_iCre1355_auto_GAPFILL.xml', help='Path to the input SBML model file.')
    # parser.add_argument('--sbml_mixo', default='./data/sink/xmls/MNL_iCre1355_mixo_GAPFILL.xml', help='Path to the input SBML model file.')
    # parser.add_argument('--sbml_hetero', default='./data/sink/xmls/MNL_iCre1355_hetero_GAPFILL.xml', help='Path to the input SBML model file.')
    # Objectives for the Trophic Modes
    parser.add_argument('--obj_auto', default='Biomass_Chlamy_auto', help='Objective reaction for growth experiment.')
    parser.add_argument('--obj_mixo', default='Biomass_Chlamy_mixo', help='Objective reaction for growth experiment.')
    parser.add_argument('--obj_hetero', default='Biomass_Chlamy_hetero', help='Objective reaction for growth experiment.')
    # Output Directory
    parser.add_argument('-o', '--outdir', default="./res/growth", help='Directory to save the growth experiment results.')
    args = parser.parse_args()

    # Model paths based on model type
    args.sbml_auto = f'./data/{args.model_type}/xmls/iCre1355_auto.xml'
    args.sbml_mixo = f'./data/{args.model_type}/xmls/iCre1355_mixo.xml'
    args.sbml_hetero = f'./data/{args.model_type}/xmls/iCre1355_hetero.xml'
    args.outdir = os.path.join(args.outdir, args.model_type)

    # Create Directory
    os.makedirs(args.outdir, exist_ok=True)

    # Debug
    print("Model import from SBML files... \n{}, \n{}, \n{}".format(args.sbml_auto, args.sbml_mixo, args.sbml_hetero))

    # Model Import
    try:
        ref_auto = io.read_sbml_model(args.sbml_auto)
    except Exception as e:
        print(f"Error loading auto model: {e}")
        sys.exit(1)

    try:
        ref_mixo = io.read_sbml_model(args.sbml_mixo)
    except Exception as e:
        print(f"Error loading mixo model: {e}")
        sys.exit(1)

    try:
        ref_hetero = io.read_sbml_model(args.sbml_hetero)
    except Exception as e:
        print(f"Error loading hetero model: {e}")
        sys.exit(1)

    # print(type(ref_auto), type(ref_mixo), type(ref_hetero))
    # print(ref_auto.objective, ref_mixo.objective, ref_hetero.objective)
    # print(ref_auto.slim_optimize(), ref_mixo.slim_optimize(), ref_hetero.slim_optimize())

    # for rxn_id in ['EX_co2_e', 'EX_nh4_e', 'EX_pi_e', 'EX_ac_e']:
    #     r = ref_auto.reactions.get_by_id(rxn_id)
    #     print(f"{rxn_id}: default bounds = {r.bounds}")
    # exit(1)

    # Run growth experiment
    experiments = {
        'Autotrophic_Rep1': GrowthExperimentItem(
            model=ref_auto,
            mu_obs=0.033,
            objective='Biomass_Chlamy_auto',
            co2_ub=1.289,
            nh4_ub=-0.569,
            pi_ub=-0.032,
            acet_ub='NA', 
            light_ub=None,
        ),
        'Autotrophic_Rep2': GrowthExperimentItem(
            model=ref_auto,
            mu_obs=0.033,
            objective='Biomass_Chlamy_auto',
            co2_ub='ND',
            nh4_ub=-0.734,
            pi_ub=-0.024,
            acet_ub='NA', 
            light_ub=None,
        ),
        'Autotrophic_Rep3': GrowthExperimentItem(
            model=ref_auto,
            mu_obs=0.042,
            objective='Biomass_Chlamy_auto',
            co2_ub=1.752,
            nh4_ub=-0.514,
            pi_ub=-0.025,
            acet_ub='NA',
            light_ub=None,
        ),
        'Autotrophic_Rep4': GrowthExperimentItem(
            model=ref_auto,
            mu_obs=0.042,
            objective='Biomass_Chlamy_auto',
            co2_ub=1.775,
            nh4_ub=-0.486,
            pi_ub=-0.026,
            acet_ub='NA',
            light_ub=None,
        ),
        'Mixotrophic_Rep1': GrowthExperimentItem(
            model=ref_mixo,
            mu_obs=0.082,
            objective='Biomass_Chlamy_mixo',
            co2_ub='ND',
            nh4_ub=-3.035,
            pi_ub=-0.099,
            acet_ub=-1.364,
            light_ub=None,
        ),
        'Mixotrophic_Rep2': GrowthExperimentItem(
            model=ref_mixo,
            mu_obs=0.082,
            objective='Biomass_Chlamy_mixo',
            co2_ub='ND',
            nh4_ub=-0.748,
            pi_ub=-0.078,
            acet_ub=-2.261,
            light_ub=None,
        ),
        'Mixotrophic_Rep3': GrowthExperimentItem(
            model=ref_mixo,
            mu_obs=0.061,
            objective='Biomass_Chlamy_mixo',
            co2_ub='ND',
            nh4_ub=-1.259,
            pi_ub=-0.047,
            acet_ub=-1.059,
            light_ub=None,
        ),
        'Mixotrophic_Rep4': GrowthExperimentItem(
            model=ref_mixo,
            mu_obs=0.061,
            objective='Biomass_Chlamy_mixo',
            co2_ub='ND',
            nh4_ub=-1.071,
            pi_ub=-0.041,
            acet_ub=-1.234,
            light_ub=None,
        ),
        'Heterotrophic_Rep1': GrowthExperimentItem(
            model=ref_hetero,
            mu_obs=0.024,
            objective='Biomass_Chlamy_hetero',
            co2_ub='ND',
            nh4_ub=-0.439,
            pi_ub=-0.010,
            acet_ub=-1.579,
            light_ub=0.0,
        ),
        'Heterotrophic_Rep2': GrowthExperimentItem(
            model=ref_hetero,
            mu_obs=0.028,
            objective='Biomass_Chlamy_hetero',
            co2_ub='ND',
            nh4_ub=-0.454,
            pi_ub=-0.012,
            acet_ub=-1.599,
            light_ub=0.0,
        ),
    }
    growth_experiment(experiments, args.outdir)
