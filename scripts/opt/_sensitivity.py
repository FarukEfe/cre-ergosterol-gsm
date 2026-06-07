from datetime import datetime
from pyexpat import model
import os, sys, argparse, hashlib, json, pandas as pd

import cobra
from cobra import io
from cobra import Model
from scripts.opt._fba import run_flux_balance_analysis
from scripts.opt._fva import run_flux_variability_analysis

from scripts.opt.plot._plot_sweep import _plot_sweep
from scripts.opt.plot._plot_dual_variables import _plot_dual_variables

def run_wildtype_baseline(model: Model, objective: str, outdir: str) -> float:
    """
    Step 1 — Wildtype baseline FBA under biomass objective.
    Returns wt_mu for normalization in downstream steps.
    """

    # Verify biomass reaction exists
    if objective not in [r.id for r in model.reactions]:
        raise ValueError(f"Biomass reaction '{objective}' not found in model.")

    with model:
        model.objective = objective
        wt_solution = model.optimize(raise_error=False)
    
    if wt_solution is None or wt_solution.status != 'optimal':
        raise ValueError("Wildtype optimization failed. Check model feasibility and objective function.")
    
    wt_mu = wt_solution.objective_value
    print(f"Wildtype baseline growth rate (wt_mu): {wt_mu:.4f}")

    # Save full flux distribution
    fluxes = pd.Series(wt_solution.fluxes).reset_index()
    fluxes.columns = ['reaction_id', 'flux_mmol_gDW_h']
    fluxes = fluxes.sort_values('flux_mmol_gDW_h', key=abs, ascending=False)
    fluxes.to_csv(os.path.join(outdir, 'wt_baseline.csv'), index=False)

    # Sterol fluxes
    sterol_rxns = ['ERG','ERG4_7ENOL','ERG4_TERMINAL']
    print("\nSterol sink fluxes at WT baseline:")
    for sink_id in sterol_rxns:
        if sink_id in wt_solution.fluxes:
            print(f"  {sink_id}: {wt_solution.fluxes[sink_id]:.6f} mmol/gDW/h")

    # Update summary.json with wt_mu
    summary_path = os.path.join(outdir, 'summary.json')
    with open(summary_path, 'r') as f:
        meta = json.load(f)
    meta['wt_mu'] = wt_mu
    with open(summary_path, 'w') as f:
        json.dump(meta, f, indent=2)

    print(f"\nWT baseline saved to {outdir}/wt_baseline.csv")
    return wt_mu

def nutrient_uptake_sweep(
        model: Model, 
        wt_mu: float,
        target_rxns: dict[str, str],
        objective: str,
        outdir: str,
        steps: int = 30
) -> dict[str, pd.DataFrame]:
    
    all_results = {}

    for rxn_id, display_name in target_rxns.items():

        try:
            rxn = model.reactions.get_by_id(rxn_id)
        except KeyError:
            print(f"Warning: Reaction '{rxn_id}' not found in model. Skipping.")
            continue

        default_ub = rxn.upper_bound

        with model:
            fva = run_flux_variability_analysis(
                model, 
                objective=objective,
                reactions=[rxn_id]
            )
            sweep_min = 0.0
            sweep_max = float(fva.loc[rxn_id, 'maximum'])
            
        if sweep_max - sweep_min < model.tolerance:
            print(f"  WARNING: {rxn_id} has no feasible range — skipping.")
            continue

        gradient = [
            sweep_min + i * (sweep_max - sweep_min) / (steps - 1)
            for i in range(steps)
        ]

        records = []
        model.objective = objective
        for step in gradient:
            with model:
                model.reactions.get_by_id(rxn_id).upper_bound = step
                sol = model.optimize(raise_error=False)

                if sol is None or sol.status != 'optimal':
                    records.append({
                        'reaction_id':    rxn_id,
                        'nutrient':       display_name,
                        'step_ub':        step,
                        'mu':             None,
                        'normalized_mu':  None,
                        'status':         getattr(sol, 'status', 'infeasible'),
                    })
                else:
                    records.append({
                        'reaction_id':    rxn_id,
                        'nutrient':       display_name,
                        'step_ub':        step,
                        'mu':             sol.objective_value,
                        'normalized_mu':  sol.objective_value / wt_mu if wt_mu > 0 else None,
                        'status':         'optimal',
                    })
        
        df = pd.DataFrame(records)
        csv_path = os.path.join(outdir, f'sweep_{rxn_id}.csv')
        df.to_csv(csv_path, index=False)
        print(f"    Saved: {csv_path}")

        # Plot results
        _plot_sweep(
            df, rxn_id, display_name, wt_mu, outdir, 
            ylabel=None if 'Biomass_Chlamy' in objective else model.reactions.get_by_id(objective).id + " flux (µ / µ_WT)"
        )

        all_results[rxn_id] = df

    return all_results

def dual_variable_analysis(
    model: Model,
    objective: str,
    target_mets: list[str],
    target_rxns: list[str],
    pathway_groups: dict[str, list[str]],
    group_colors: dict[str, str],
    outdir: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Run FBA to get WT baseline, and shadow prices + reduced costs for relevant reactions.
    """
    # Solve FBA
    with model:
        model.objective = objective
        soln = model.optimize(raise_error=False)

    # Verify solution
    if soln is None or soln.status != 'optimal':
        raise ValueError("Optimization failed. Check model feasibility and objective function.")
    
    # Extract duals
    sp_df = soln.shadow_prices.dropna()
    sp_df = sp_df[sp_df.index.isin(target_mets)]
    rc_df = soln.reduced_costs.dropna()
    rc_df = rc_df[rc_df.index.isin(target_rxns)]

    sp_df = sp_df.reset_index()
    sp_df.columns = ['metabolite_id', 'shadow_price']
    rc_df = rc_df.reset_index()
    rc_df.columns = ['reaction_id', 'reduced_cost']

    # Rank by bottleneck
    sp_df['abs_shadow_price'] = sp_df['shadow_price'].abs()
    sp_df = sp_df.sort_values('abs_shadow_price', ascending=False).reset_index(drop=True)
    sp_df['rank'] = sp_df.index + 1

    rc_df['abs_reduced_cost'] = rc_df['reduced_cost'].abs()
    rc_df = rc_df.sort_values('abs_reduced_cost', ascending=False).reset_index(drop=True)
    rc_df['rank'] = rc_df.index + 1

    # Save results
    sp_df.to_csv(os.path.join(outdir, 'shadow_prices.csv'), index=False)
    rc_df.to_csv(os.path.join(outdir, 'reduced_costs.csv'), index=False)

    # Update summary.json
    with open(os.path.join(outdir, 'summary.json'), 'r') as f:
        meta = json.load(f)
    meta['top_shadow_price_mets'] = sp_df.head(10)['metabolite_id'].tolist()
    meta['top_reduced_cost_rxns'] = rc_df.head(10)['reaction_id'].tolist()
    with open(os.path.join(outdir, 'summary.json'), 'w') as f:
        json.dump(meta, f, indent=2)
    
    # Plot dual variable results
    _plot_dual_variables(sp_df, rc_df, pathway_groups, group_colors, outdir)

    return sp_df, rc_df



if __name__ == "__main__":

    # Script Argument(s)
    parser = argparse.ArgumentParser(
        prog='sensitivity',
        description='Run sensitivity analysis on your model by performing FBA and FVA with different objective functions.'
    )
    parser.add_argument('--sbmlpath', default="./data/sink/xmls/MNL_iCre1355_auto_GAPFILL.xml", help='Path to the input SBML model file.')
    parser.add_argument('--objective', default="Biomass_Chlamy_auto", help='Objective function to optimize for sensitivity analysis.')
    parser.add_argument('--outdir', default="./res/sensitivity", help='Directory to save the sensitivity analysis results.')
    parser.add_argument('--steps', type=int, default=30, help='Number of steps for nutrient uptake sweep.')
    args = parser.parse_args()

    args.outdir = os.path.join(args.outdir, os.path.splitext(os.path.basename(args.sbmlpath))[0], args.objective)
    os.makedirs(args.outdir, exist_ok=True)

    # Load model
    ref, err = io.validate_sbml_model(args.sbmlpath)
    if not ref:
        print(f"Error loading model: {err}")
        sys.exit(1)

    # Slim optimize
    ref.objective = args.objective
    ref.tolerance = 1e-6
    ref.slim_optimize()

    # Log hash and cobra version
    model_md5 = hashlib.md5(open(args.sbmlpath,'rb').read()).hexdigest()
    timestamp = datetime.now().isoformat()

    # Create metadata
    meta = {
        'timestamp':       timestamp,
        'cobраpy_version': cobra.__version__,
        'solver':          ref.solver.name,
        'model_md5':       model_md5,
        'sbmlpath':        args.sbmlpath,
        'n_reactions':     len(ref.reactions),
        'n_metabolites':   len(ref.metabolites),
        'n_genes':         len(ref.genes),
    }

    with open(os.path.join(args.outdir, 'summary.json'), 'w') as f:
        json.dump(meta, f, indent=2)

    # Wildtype baseline
    wt_mu = run_wildtype_baseline(model=ref, objective=args.objective, outdir=args.outdir)

    # Nutrient sweep sensitivity analysis
    print("Running nutrient uptake sweep ...\n")
    uptake_rxns = {
        'PRISM_design_growth': 'Photon',
        'CO2t':                'CO2',
        'NH4t':                'Ammonium',
        'O2t':                 'Oxygen',
        'PIt':                 'Phosphate',
        'FE2GTPabc':           'Iron(II)',
        'FE3t':                'Iron(III)',
    }

    nutrient_uptake_results = nutrient_uptake_sweep(
        model=ref,
        wt_mu=wt_mu,
        target_rxns=uptake_rxns,
        objective=args.objective,
        outdir=args.outdir,
        steps=args.steps
    )

    # LP dual sensitivity analysis
    print("\nRunning dual variable analysis ...\n")

    # METABOLITES

    mep_pathway_metabolites = [
        # MEP pathway intermediates (chloroplast)
        'pyr_h',      # pyruvate (chloroplast)
        'g3p_h',      # glyceraldehyde-3-phosphate (chloroplast)
        'dxyl5p_h',   # 1-deoxy-D-xylulose-5-phosphate (DXP)
        'mep_h',      # 2-C-methyl-D-erythritol-4-phosphate (MEP)
        'cdpme_h',    # CDP-ME
        'cdpmep_h',   # CDP-MEP
        'mecpp_h',    # MEcPP
        'hmbpp_h',    # HMBPP
        'ipdp_h',     # IPP (chloroplast)
        'dmapp_h',    # DMAPP (chloroplast)
    ]

    isoprenoid_assembly_metabolites = [
        'ipdp_c',     # IPP (cytosol)
        'dmapp_c',    # DMAPP (cytosol)
        'grdp_c',     # GPP
        'frdp_c',     # FPP
        'ggdp_c',     # GGPP
        'sqlne_c',    # squalene
        'sqlnol_c',   # 2,3-oxidosqualene
    ]

    sterol_trunk_metabolites = [
        'cyartenol_c',
        'obfool_c',
        '4amethfec_c',
        '24menlophenol_c',
        'episterol_c',
        'ergtrienol_c',
    ]

    sterol_branch_metabolites = [
        'ergosta57dienol_c',
        'ergosta5722trienol_c',
        'ergosterol_c',
        'ergost7enol_c',
        '7dhporiferasterol_c',
    ]

    mva_pathway_metabolites = [
        'hmgcoa_c',
        'mva_c',
        'mva5p_c',
        'diphmva_c',
    ]

    medium_metabolites = [
        'photonVis_e',
        'co2_e',
        'nh4_e',
        'o2_e',
        'pi_e',
        'fe2_e',
        'fe3_e',
    ]

    # Combined for filtering
    relevant_metabolites = (
        mep_pathway_metabolites +
        isoprenoid_assembly_metabolites +
        sterol_trunk_metabolites +
        sterol_branch_metabolites +
        mva_pathway_metabolites +
        medium_metabolites
    )

    pathway_groups = {
        "Sterol branch": sterol_branch_metabolites,
        "Sterol trunk": sterol_trunk_metabolites,
        "Isoprenoid / MEP": isoprenoid_assembly_metabolites + mep_pathway_metabolites,
        "MVA pathway": mva_pathway_metabolites,
        "Medium": medium_metabolites,
    }

    # REACTIONS

    # Trunk (shared by both branches)
    trunk_reactions = [
        'SMO',           # squalene → 2,3-oxidosqualene
        'CAS',           # 2,3-oxidosqualene → cycloartenol
        'CYARTMT',       # cycloartenol → 24-methylenecycloartenol
        'CYEUOLS',       # → cycloeucalenol
        'CYEUOLCYCI',    # cycloeucalenol → obtusifoliol
        'OBFOOLOR',      # obtusifoliol → mergtrol
        'MERGTROLR',     # mergtrol → 4α-methylfecosterol
    ]

    # Pre-branch isomerization
    prebranch_reactions = [
        'HYD1',          # mfecostrl → methylop_c (ergosterol branch entry)
        'CDI1_BRANCH',   # mfecostrl → 24menlophenol_c (7DHP branch entry)
    ]

    # Ergosterol branch
    ergosterol_reactions = [
        'ERG28',         # methylop_c → episterol
        'ERG3',          # episterol → ergtrienol
        'ERG5',          # ergtrienol → ergtetraenol
        'ERG4_ERGSTRL',           # ergtetraenol → ergosterol
        'ERG4_7ENOL',    # ergtrienol → ergost7enol (minor dead-end)
    ]

    # 7-DHPoriferasterol branch
    dhp_reactions = [
        'SMT1_2ND',      # 24menlophenol → 24ethlophenol
        'C4DEMETH_BRANCH',# 24ethlophenol → avenosterol
        'SMT1_CROSS',    # episterol → avenosterol (shortcut)
        'STE1_ERG3',     # avenosterol → stigmatrienol
        'CYP710_ERG5',   # stigmatrienol → stigmatrienol22
        'ERG4_7DEHYD', # stigmatrienol22 → 7dhporiferasterol
    ]

    # Combined
    sterol_pathway_reactions = (
        trunk_reactions + 
        prebranch_reactions + 
        ergosterol_reactions + 
        dhp_reactions
    )

    dual_variable_analysis(
        model=ref,
        objective=args.objective,
        target_mets=relevant_metabolites,
        target_rxns=sterol_pathway_reactions,
        pathway_groups=pathway_groups,
        group_colors={
            "Sterol branch":    "#d95f02",
            "Sterol trunk":     "#e6ab02",
            "Isoprenoid / MEP": "#1b9e77",
            "MVA pathway":      "#7570b3",
            "Medium":           "#66a61e",
            "Other":            "#aaaaaa",
        },
        outdir=args.outdir
    )

    print("\nSensitivity analysis complete.")






