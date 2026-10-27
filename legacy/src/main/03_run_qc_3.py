from pathlib import Path
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## create the dataframe
tinception_dir = Path("/Volumes/Extreme_SSD/payam_data/Tinception")
qc_dir = tinception_dir / "qc_dir"

my_dict = {
            "tinception_id": [],
            "cjv": [],
            "cnr": [],
            'snr_csf': [],
            'snr_gm': [],
            'snr_wm': [],
            'snr_total': []
            }


for fname_dir in tqdm(list(qc_dir.iterdir())):
    name = fname_dir.name
    if name.startswith(".") or name == "logs" or not fname_dir.is_dir():
        continue
    else:
        df_json = pd.read_json(fname_dir / "anat" / f"{name}_T1w.json")
        my_dict["tinception_id"].append(name)
        for key in list(my_dict.keys())[1:]:
            my_dict[key].append(df_json[key].unique()[0])


df = pd.DataFrame(my_dict)
df = pd.melt(
            df,
            id_vars=['tinception_id'],
            value_vars=['cjv', 'cnr', 'snr_csf', 'snr_gm', 'snr_total', 'snr_wm'],
            var_name='method',
            value_name='value'
            )

df_master = pd.read_csv("../../master_files/master.csv")
df = df.merge(
                df_master[["tinception_id", "group", "sex", "age", "site"]],
                on="tinception_id",
                how="inner"
)
df.sort_values(by=["tinception_id", "method"], inplace=True)
df.reset_index(drop=True, inplace=True)

## plot histograms
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), "axes.linewidth": 1.5})
plt.rcParams['axes.linewidth'] = 2
custom_palette = {"CO": "#4575b4", "TI": "#d73027"} 
col_order = ["cjv", "cnr", "snr_total", "snr_gm", "snr_wm", "snr_csf"]

g = sns.FacetGrid(
    data=df,
    col="method",
    col_wrap=3,
    col_order=col_order,
    hue="group",
    palette=custom_palette,
    sharex=False,
    sharey=True,
    height=4,
    aspect=1.2,
    despine=True
)

g.map_dataframe(
    sns.histplot,
    x="value",
    fill=False,
    linewidth=1,
    kde=True,
    element="bars", 
    alpha=1,
    bins=45,
    common_norm=False,
    kde_kws={"bw_adjust": 0.8, "cut": 3}
)

for i, ax in enumerate(g.axes.flat):
    if i not in [0, 3]:
        ax.spines['left'].set_visible(False)
g.set_titles(col_template="{col_name}", weight='bold', size=14)
g.set_axis_labels("Value", "Density", fontsize=12)

g.add_legend(title="Study Group", adjust_subtitles=True, loc="upper right")
plt.tight_layout()
g.figure.savefig(
                tinception_dir / "plots" / "QC.pdf", 
                format="pdf",
                dpi=300,       
                bbox_inches="tight"
                )


## Multivariate Desicion
df_wide = df.pivot(index="tinception_id", columns="method", values="value").reset_index()

def find_outliers(data, col, direction='both'):
    q1 = data[col].quantile(0.25)
    q3 = data[col].quantile(0.75)
    iqr = q3 - q1
    low, high = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    
    if direction == 'low': return data['tinception_id'][data[col] < low]
    if direction == 'high': return data['tinception_id'][data[col] > high]
    return data['tinception_id'][(data[col] < low) | (data[col] > high)]

# Define quality directions
outliers = {
    "cjv": find_outliers(df_wide, "cjv", "high"),
    "cnr": find_outliers(df_wide, "cnr", "low"),
    "snr_total": find_outliers(df_wide, "snr_total", "low")
}

# Combine and count mentions
from collections import Counter
all_bad_ids = [idx for sublist in outliers.values() for idx in sublist]
drop_counts = Counter(all_bad_ids)

# Subjects appearing in 2+ 'bad' categories
to_drop = [ids for ids, count in drop_counts.items() if count >= 2]
print("******** Multivariate Decision ***********")
print(f"Recommended for exclusion (Multi-metric outliers): {to_drop}\n")

## Univariate Decision
outlier_summary = []

for method in col_order:
    m_data = df[df['method'] == method]
    q1, q3 = m_data['value'].quantile(0.25), m_data['value'].quantile(0.75)
    iqr = q3 - q1
    
    # Define bounds based on metric type
    if method == 'cjv':
        cutoff = q3 + 1.5 * iqr
        outliers = m_data[m_data['value'] > cutoff]
    else:
        cutoff = q1 - 1.5 * iqr
        outliers = m_data[m_data['value'] < cutoff]
        
    for _, row in outliers.iterrows():
        outlier_summary.append({
            "ID": row['tinception_id'], 
            "Metric": method, 
            "Value": round(row['value'], 3),
            "Site": row['site'],
            "Group": row['group']
        })

outlier_report = pd.DataFrame(outlier_summary)
print("******** Univariate Decision ***********")
print(outlier_report.sort_values(by="ID").reset_index(drop=True))
