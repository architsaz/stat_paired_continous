import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
import shutil
import sys

dir_runfebio = "/dagon1/achitsaz/runfebio/" 
table_dir = "./table/"
save_dir = "./figs/"
list_of_region = ["aneu", "dome", "body", "neck", "thin", "thick", "rupt", "bleb"]
    # clean save data directory
for region in list_of_region:
    region_figs = f'{save_dir}{region}'
    if os.path.exists(region_figs):
        shutil.rmtree(region_figs)
    os.makedirs(region_figs, exist_ok=True) 

dfs = {
    "aneu" : pd.read_csv(f"{table_dir}aneu/continous_fields.csv"),
    "dome" : pd.read_csv(f"{table_dir}dome/continous_fields.csv"),
    "body" : pd.read_csv(f"{table_dir}body/continous_fields.csv"),
    "neck" : pd.read_csv(f"{table_dir}neck/continous_fields.csv"),
    "thin" : pd.read_csv(f"{table_dir}thin/continous_fields.csv"),
    "thick" : pd.read_csv(f"{table_dir}thick/continous_fields.csv"),
    "rupt" : pd.read_csv(f"{table_dir}rupt/continous_fields.csv"),
    "bleb" : pd.read_csv(f"{table_dir}bleb/continous_fields.csv")
}
# Plotting bar figure for each continous fiels with different statisctic metrics
for region in list_of_region:
    save_path = save_dir+region
    # Labels and values
    lable = ["mean_von_homo", "mean_von_hete",
            "max_von_homo", "max_von_hete",
            "mean_UQ_von_homo", "mean_UQ_von_hete"
            ]
    means = [dfs[region]["mean_von_homo"].mean(), dfs[region]["mean_von_hete"].mean(),
            dfs[region]["max_von_homo"].mean(), dfs[region]["max_von_hete"].mean(),
            dfs[region]["mean_UQ_von_homo"].mean(), dfs[region]["mean_UQ_von_hete"].mean(),
            ]
    stds = [dfs[region]["mean_von_homo"].std(), dfs[region]["mean_von_hete"].std(),
            dfs[region]["max_von_homo"].std(), dfs[region]["max_von_hete"].std(),
            dfs[region]["mean_UQ_von_homo"].std(), dfs[region]["mean_UQ_von_hete"].std(),
            ]

    # Create bar plot with error bars
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.bar(lable, means, yerr=stds, capsize=7, color=['#ff9999', '#66b3ff','#ff9999', '#66b3ff','#ff9999', '#66b3ff'])
    ax.set_xticklabels(lable, rotation=20)
    # Labels and title
    ax.set_title("Comparison Von Mises Bewtween Hete. and Homo. Conditions")
    ax.set_ylabel("Von Mises Stress (dyn/cm²)")
    # Custom legend
    legend_elements = [
        Patch(facecolor='#ff9999', label='Homogeneous'),
        Patch(facecolor='#66b3ff', label='Heterogeneous')
    ]
    ax.legend(handles=legend_elements)
    # Save figure
    plot_file = os.path.join(save_path, "von_plot.png")
    print(plot_file)
    plt.savefig(plot_file)
    plt.close()

    # Labels and values
    lable = ["mean_eval_max_homo", "mean_eval_max_hete",
            "max_eval_max_homo", "max_eval_max_hete",
            "mean_UQ_eval_max_homo", "mean_UQ_eval_max_hete"
            ]
    means = [dfs[region]["mean_eval_max_homo"].mean(), dfs[region]["mean_eval_max_hete"].mean(),
            dfs[region]["max_eval_max_homo"].mean(), dfs[region]["max_eval_max_hete"].mean(),
            dfs[region]["mean_UQ_eval_max_homo"].mean(), dfs[region]["mean_UQ_eval_max_hete"].mean(),
            ]
    stds = [dfs[region]["mean_eval_max_homo"].std(), dfs[region]["mean_eval_max_hete"].std(),
            dfs[region]["max_eval_max_homo"].std(), dfs[region]["max_eval_max_hete"].std(),
            dfs[region]["mean_UQ_eval_max_homo"].std(), dfs[region]["mean_UQ_eval_max_hete"].std(),
            ]

    # Create bar plot with error bars
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.bar(lable, means, yerr=stds, capsize=7, color=['#ff9999', '#66b3ff','#ff9999', '#66b3ff','#ff9999', '#66b3ff'])
    ax.set_xticklabels(lable, rotation=20)
    # Labels and title
    ax.set_title("Comparison Max of Eigenvalue Bewtween Hete. and Homo. Conditions")
    ax.set_ylabel("Max of Eigenvalue (dyn/cm²)")
    # Custom legend
    legend_elements = [
        Patch(facecolor='#ff9999', label='Homogeneous'),
        Patch(facecolor='#66b3ff', label='Heterogeneous')
    ]
    ax.legend(handles=legend_elements)
    # Save figure
    plot_file = os.path.join(save_path, "eval_max_plot.png")
    print(plot_file)
    plt.savefig(plot_file)
    plt.close()

    # Labels and values
    lable = ["concen_von_homo", "concen_von_hete",
            "concen_eval_max_homo", "concen_eval_max_hete",
            "force_ratio_homo", "force_ratio_hete",
            "area_ratio_homo", "area_ratio_hete"
            ]
    means = [dfs[region]["concen_von_homo"].mean(), dfs[region]["concen_von_hete"].mean(),
            dfs[region]["concen_eval_max_homo"].mean(), dfs[region]["concen_eval_max_hete"].mean(),
            dfs[region]["force_ratio_homo"].mean(), dfs[region]["force_ratio_hete"].mean(),
            dfs[region]["area_ratio_homo"].mean(), dfs[region]["area_ratio_hete"].mean()
            ]
    stds = [dfs[region]["concen_von_homo"].std(), dfs[region]["concen_von_hete"].std(),
            dfs[region]["concen_eval_max_homo"].std(), dfs[region]["concen_eval_max_hete"].std(),
            dfs[region]["force_ratio_homo"].std(), dfs[region]["force_ratio_hete"].std(),
            dfs[region]["area_ratio_homo"].std(), dfs[region]["area_ratio_hete"].std()        
            ]

    # Create bar plot with error bars
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.bar(lable, means, yerr=stds, capsize=7, color=['#ff9999', '#66b3ff','#ff9999', '#66b3ff','#ff9999', '#66b3ff','#ff9999', '#66b3ff'])
    ax.set_xticklabels(lable, rotation=20)
    # Labels and title
    ax.set_title("Comparison Concentration of Stress Bewtween Hete. and Homo. Conditions")
    ax.set_ylabel("Ratio")
    # Custom legend
    legend_elements = [
        Patch(facecolor='#ff9999', label='Homogeneous'),
        Patch(facecolor='#66b3ff', label='Heterogeneous')
    ]
    ax.legend(handles=legend_elements)
    # Save figure
    plot_file = os.path.join(save_path, "concen_stress_plot.png")
    print(plot_file)
    plt.savefig(plot_file)
    plt.close()