import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def _plot_dual_variables(sp_df, rc_df, pathway_groups, group_colors, outdir):
    
    _POS_COL  = "#2c7bb6"
    _NEG_COL  = "#d7191c"
    _ZERO_COL = "#aaaaaa"
    
    def _met_group(met_id):
        for group, members in pathway_groups.items():
            if met_id in members:
                return group
        return "Other"


    def _sign_colour(val):
        if val > 0:
            return _POS_COL
        elif val < 0:
            return _NEG_COL
        return _ZERO_COL
    
    TOP_N      = 20
    NOISE_SP   = 1e-6
    NOISE_RC   = 1e-9

    # Shadow Prices
    sp = sp_df.dropna(subset=["shadow_price"]).copy()
    sp["abs_shadow_price"] = sp["shadow_price"].abs()
    sp = sp[sp["abs_shadow_price"] > NOISE_SP]
    sp = sp.sort_values("abs_shadow_price", ascending=False).head(TOP_N)
    sp = sp.sort_values("shadow_price", ascending=True)
    sp["group"] = sp["metabolite_id"].apply(_met_group)

    # Reduced Costs
    rc = rc_df.dropna(subset=["reduced_cost"]).copy()
    rc["abs_reduced_cost"] = rc["reduced_cost"].abs()
    rc = rc[rc["abs_reduced_cost"] > NOISE_RC]
    rc = rc.sort_values("abs_reduced_cost", ascending=False).head(TOP_N)
    rc = rc.sort_values("reduced_cost", ascending=True)

    # Figure Setup
    fig, (ax_sp, ax_rc) = plt.subplots(
        1, 2,
        figsize=(16, max(6, TOP_N * 0.38)),
        constrained_layout=True,
    )
    fig.patch.set_facecolor("white")

    # Panel A — shadow prices (coloured by pathway group)
    ax_sp.set_facecolor("white")
    colours = [group_colors[g] for g in sp["group"]]
    ax_sp.barh(sp["metabolite_id"], sp["shadow_price"],
               color=colours, edgecolor="white", linewidth=0.4, height=0.72)
    ax_sp.axvline(0, color="#333333", linewidth=0.8, linestyle="--")
    ax_sp.set_xlabel("Shadow price  (Δ obj / Δ metabolite balance)", fontsize=9)
    ax_sp.set_title("A  Metabolite shadow prices", fontsize=11, fontweight="bold", loc="left")
    ax_sp.tick_params(axis="y", labelsize=7.5)
    ax_sp.tick_params(axis="x", labelsize=8)
    ax_sp.spines[["top", "right"]].set_visible(False)

    present = sp["group"].unique()
    legend_patches = [
        mpatches.Patch(facecolor=group_colors[g], label=g)
        for g in group_colors if g in present
    ]
    if legend_patches:
        ax_sp.legend(handles=legend_patches, fontsize=7, loc="lower right",
                     framealpha=0.7, title="Pathway", title_fontsize=7)

    # Panel B — reduced costs (coloured by sign)
    ax_rc.set_facecolor("white")
    rc_colours = [_sign_colour(v) for v in rc["reduced_cost"]]
    ax_rc.barh(rc["reaction_id"], rc["reduced_cost"],
               color=rc_colours, edgecolor="white", linewidth=0.4, height=0.72)
    ax_rc.axvline(0, color="#333333", linewidth=0.8, linestyle="--")
    ax_rc.set_xlabel("Reduced cost  (Δ obj / forcing reaction flux)", fontsize=9)
    ax_rc.set_title("B  Reaction reduced costs", fontsize=11, fontweight="bold", loc="left")
    ax_rc.tick_params(axis="y", labelsize=7.5)
    ax_rc.tick_params(axis="x", labelsize=8)
    ax_rc.spines[["top", "right"]].set_visible(False)

    sign_patches = [
        mpatches.Patch(facecolor=_POS_COL, label="Positive (beneficial to force)"),
        mpatches.Patch(facecolor=_NEG_COL, label="Negative (costly to force)"),
    ]
    ax_rc.legend(handles=sign_patches, fontsize=7, loc="lower right",
                 framealpha=0.7, title="Direction", title_fontsize=7)
    
    # Save figure
    out_path = os.path.join(outdir, "dual_variables.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out_path}")