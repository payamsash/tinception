from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap

from brainspace.gradient import GradientMaps
from brainspace.utils.parcellation import map_to_labels
from brainspace.datasets import load_parcellation
from neuromaps.datasets import fetch_fslr
from surfplot import Plot


## paths
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
plots_dir = tinception_dir / "plots" / "ssa"
ssa_results_dir = tinception_dir / "ssa_results"
df_lambda = pd.read_csv(ssa_results_dir / "lambdas.csv")
df_grad = pd.read_csv(ssa_results_dir / "grads.csv")
plotting_kwargs = {
                "format": "pdf",
                "dpi": 300,
                "bbox_inches": "tight"
}

## plotting scree plot
palette_color = ["#145e92", "#b51f1f"]
g = sns.FacetGrid(data=df_lambda, aspect=2.5, height=3, ylim=[0, 18])
g.map_dataframe(
    sns.stripplot,
    x='component',
    y='value',
    hue="group",
    dodge=True,
    palette=palette_color,
    hue_order=["CO", "TI"],
    size=2.5,
    alpha=0.15
    )
g.map_dataframe(
    sns.boxplot,
    x='component',
    y='value',
    hue="group",
    dodge=True,
    width=0.8,
    gap=0.2,
    color="white",
    fill=False,
    showfliers=False,
    palette=palette_color,
    hue_order=["CO", "TI"]
    )

g.refline(x=2.5)
g.add_legend(frameon=False, bbox_to_anchor=(0.9, 0.9))
g.savefig(
            plots_dir / "lambda.pdf",
            **plotting_kwargs
        )

## plot state space plots of average conns
conn_co = np.load(ssa_results_dir / "avg_connectome_co.npy", allow_pickle=False)
conn_ti = np.load(ssa_results_dir / "avg_connectome_ti.npy", allow_pickle=False)

gm = GradientMaps(
                    n_components=3,
                    approach="dm",
                    kernel=None,
                    alignment='procrustes',
                    random_state=0
                    )
gm.fit([conn_co, conn_ti])

gm_vals = np.array(gm.aligned_).flatten()
group_labels = np.repeat(['CO', 'TI'], 1000 * 3)
component_labels = np.tile([1, 2, 3], 2 * 1000)
df_ss = pd.DataFrame({
                'group': group_labels,
                'component': component_labels,
                'value': gm_vals
            })

fname_lut = ssa_results_dir / "Schaefer2018_1000Parcels_7Networks_order.txt"
df_lut = pd.read_csv(fname_lut, delimiter="\t", header=None, names=["x", "labels", "z", "a", "b", "c"])
df_lut = df_lut[["x", "labels"]]

network_map = {
    "Vis": "Visual",
    "SomMot": "Somatomotor",
    "DorsAttn": "Dorsal Attention",
    "SalVentAttn": "Ventral Attention",
    "Limbic": "Limbic",
    "Cont": "Frontoparietal",
    "Default": "Default"
}

def get_network(label):
    for key, value in network_map.items():
        if key in label:
            return value
    return "Unknown"

df_lut['Network'] = df_lut['labels'].apply(get_network)

## add networks
parcel_indices = np.repeat(np.arange(1, 1001), 3) # 3000 entries
full_parcel_column = np.tile(parcel_indices, 2)   # 6000 entries total
df_ss['x'] = full_parcel_column
df_ss = df_ss.merge(df_lut[['x', 'Network', 'labels']], on='x', how='left')
df_ss.drop(columns='x', inplace=True)
df_wide = df_ss.pivot(
    index=['group', 'labels', 'Network'], 
    columns='component', 
    values='value'
).reset_index()
df_wide.columns = [f'Comp_{c}' if isinstance(c, int) else c for c in df_wide.columns]
df_ss = df_wide.melt(
    id_vars=['group', 'labels', 'Network', 'Comp_1'], 
    value_vars=['Comp_2', 'Comp_3'], 
    var_name='other_component', 
    value_name='other_value'
)

palette = 5 * ["gainsboro"] + [
                sns.cubehelix_palette(rot=-.2).as_hex()[2],
                sns.cubehelix_palette(rot=.2).as_hex()[2]
                ]
hue_order = ['Dorsal Attention', 'Ventral Attention', 'Limbic', 'Frontoparietal', 'Visual', 'Default', 'Somatomotor']

g = sns.relplot(
    data=df_ss,
    x="Comp_1",            
    y="other_value",       
    row="other_component", 
    col="group",           
    kind="scatter",
    hue="Network",
    hue_order=hue_order,
    palette=palette,
    height=3.5,
    aspect=1.6,
    facet_kws={'sharex': True, 'sharey': False},
    color="white",
    edgecolor="k",
    alpha=0.79,
    legend=False,
)
for i in range(2):
    g.axes[i][1].spines['left'].set_visible(False)
    g.axes[i][1].tick_params(axis='y', left=False, labelleft=False)

g.savefig(
            plots_dir / "state_space.pdf",
            **plotting_kwargs
        )


## plot network schaefer
surfaces = fetch_fslr()
lh, rh = surfaces['inflated']
p = Plot(lh, rh, background=(0, 0, 1), brightness=0.65)
lh_parc, rh_parc = load_parcellation('schaefer', scale=1000)

region_numbers = df_lut['labels'].values.tolist()
regions = np.where(np.isin(lh_parc, region_numbers), lh_parc, 0)


alpha = 0.79
color_other = 'k'
cmap_other = ListedColormap(color_other, 'regions', N=1)

color_dmn = palette[-2]
cmap_dmn = ListedColormap(color_dmn, 'regions', N=1)
dmn_ids = df_lut[df_lut['Network'] == 'Default']['x'].values
dmn_lh = np.where(np.isin(lh_parc, dmn_ids), lh_parc, 0)
dmn_rh = np.where(np.isin(rh_parc, dmn_ids), rh_parc, 0)

color_som = palette[-1]
cmap_som = ListedColormap(color_som, 'regions', N=1)
som_ids = df_lut[df_lut['Network'] == 'Somatomotor']['x'].values
som_lh = np.where(np.isin(lh_parc, som_ids), lh_parc, 0)
som_rh = np.where(np.isin(rh_parc, som_ids), rh_parc, 0)

p.add_layer({'left': lh_parc, 'right': rh_parc}, cmap=cmap_other, alpha=0.5,
            as_outline=True, cbar=False)

p.add_layer({'left': dmn_lh, 'right': dmn_rh}, 
            cmap=cmap_dmn, alpha=alpha,
            cbar=False)
p.add_layer({'left': som_lh, 'right': som_rh}, 
            cmap=cmap_som, alpha=alpha,
            cbar=False)

fig = p.build()
fig.tight_layout()
fig.savefig(
            plots_dir / "brain_schaefer_networks.pdf",
            **plotting_kwargs
        )


## plot brain gradients

## compute gradients
datas = []
for group in ["co", "ti"]:
    datas.append(np.load(ssa_results_dir / f"avg_connectome_{group}.npy"))

gm = GradientMaps(
                    n_components=3,
                    approach="dm",
                    kernel=None,
                    alignment="procrustes",
                    random_state=0
                    )
gm.fit(datas)
gm_co_grad_1 = gm.aligned_[0][:, 0]
gm_ti_grad_1 = gm.aligned_[0][:, 0]
gm_diff_grad_2 = gm.aligned_[1][:, 1] - gm.aligned_[0][:, 1]
gm_diff_grad_3 = gm.aligned_[1][:, 2] - gm.aligned_[0][:, 2]

## initiate brain
lh_parc, rh_parc = load_parcellation('schaefer', scale=1000)
full_parcellation = np.concatenate([lh_parc, rh_parc])
surfaces = fetch_fslr()
lh, rh = surfaces['inflated']
n_vh = lh_parc.shape[0]

## start looping to plot
gm_data = [gm_co_grad_1, gm_ti_grad_1, gm_diff_grad_2, gm_diff_grad_3]
titles = ["co_grad_1", "ti_grad_1", "diff_grad_2", "diff_grad_3"]
color_ranges = [(-25, 25), (-25, 25), (-1, 1), (-1, 1)]

for gm_d, title, color_range in zip(gm_data, titles, color_ranges):
    grad_surface_map = map_to_labels(
                                    gm_d,
                                    full_parcellation,
                                    mask=full_parcellation > 0
                                    )
    grad_lh = grad_surface_map[:n_vh]
    grad_rh = grad_surface_map[n_vh:]
    p = Plot(lh, rh, background=(0, 0, 1), brightness=0.65)
    p.add_layer(
                {'left': grad_lh, 'right': grad_rh}, 
                cmap='coolwarm',
                alpha=0.7, 
                cbar=True, color_range=color_range
                )

    p.add_layer(
                {'left': lh_parc, 'right': rh_parc}, 
                cmap=ListedColormap(['k']), 
                as_outline=True, alpha=0.5,
                cbar=False
                )

    fig = p.build()
    fig.tight_layout()
    fig.savefig(
                plots_dir / f"brain_{title}.pdf",
                **plotting_kwargs
            )
    
