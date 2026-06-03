import os, matplotlib.pyplot as plt

def _plot_sweep(df, rxn_id, display_name, wt_mu, outdir, ylabel=None):

    ylabel = ylabel if ylabel else 'Normalized growth rate (µ / µ_WT)'
    
    optimal = df[df['status'] == 'optimal']
    infeasible = df[df['status'] != 'optimal']

    fig, ax = plt.subplots(figsize=(8, 5))

    # Main curve
    ax.plot(
        optimal['step_ub'], optimal['normalized_mu'],
        'o-', linewidth=2, markersize=5,
        color='#2166ac', markerfacecolor='#fdae61',
        label='Predicted µ (normalized)'
    )

    # Mark infeasible steps
    if not infeasible.empty:
        ax.scatter(
            infeasible['step_ub'], [0] * len(infeasible),
            marker='x', color='red', s=60, zorder=5,
            label='Infeasible'
        )

    # WT baseline
    ax.axhline(y=1.0, linestyle='--', color='gray', linewidth=1.2, label='WT baseline')

    ax.set_xlabel(f'{display_name} uptake upper bound (mmol/gDW/h)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(f'Nutrient Sensitivity: {display_name}', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.25)

    plt.tight_layout()
    png_path = os.path.join(outdir, f'sweep_{rxn_id}.png')
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {png_path}")