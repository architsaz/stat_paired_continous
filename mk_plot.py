import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
import shutil
import sys

dir_runfebio = "/dagon1/achitsaz/runfebio/" 
table_dir = "./table/"
stat_dir ="./stat/"
save_dir = "./figs/"
list_of_region = ["aneu", "dome", "body", "neck", "thin", "thick", "rupt", "bleb"]
list_of_region = ["aneu"]
# Load stat data
status = "blind_status"
group = "blind_class"
stat_path = stat_dir+f"/{status}_{group}_results.csv"
try:
        stat=pd.read_csv(stat_path)
except Exception as e:
        print(f"An error occurred: {e}")         
        exit()
# clean save data directory
for region in list_of_region:
    region_figs = f'{save_dir}{region}'
    if os.path.exists(region_figs):
        shutil.rmtree(region_figs)
    os.makedirs(region_figs, exist_ok=True) 

# Load all dataframes
dfs = {region: pd.read_csv(f"{table_dir}{region}/continous_fields.csv") for region in list_of_region}


# Plotting bar figure for each continous fiels with different statisctic metrics
for region in list_of_region:
        save_path = os.path.join(save_dir, region)
        os.makedirs(save_path, exist_ok=True)

        metrics = ["mean_von", "max_von", "mean_UQ_von"]
        label = ["mean_von_homo", "mean_von_hete", "max_von_homo", "max_von_hete", "mean_UQ_von_homo", "mean_UQ_von_hete"]
        means = [dfs[region][field].mean() for field in label]
        stds = [dfs[region][field].std() for field in label]

        fig, ax = plt.subplots(figsize=(7, 7))
        bar_colors = ['#ff9999', '#66b3ff'] * len(metrics)
        ax.bar(label, means, yerr=stds, capsize=7, color=bar_colors)
        ax.set_xticks(range(len(label)))
        ax.set_xticklabels([l.replace("_homo", " (Homo)").replace("_hete", " (Hete)") for l in label], rotation=20)
        ax.set_title(f"Comparison of Von Mises Stress in {region}")
        ax.set_ylabel("Von Mises Stress (dyn/cm²)")

        legend_elements = [
                Patch(facecolor='#ff9999', label='Homogeneous'),
                Patch(facecolor='#66b3ff', label='Heterogeneous')
        ]
        ax.legend(handles=legend_elements)

        # Add asterisks if significant
        for i, metric in enumerate(metrics):
                row = stat[(stat["region"] == region) & (stat["field"] == metric)]
                if not row.empty:
                        p_val = row["p_paired_test"].values[0]
                        p_val_log = row["logit_p"].values[0]
                        odds = row["odds_ratio"].values[0]
                        ci_lower = row["odds_ratio_CI_lower"].values[0]
                        ci_upper = row["odds_ratio_CI_upper"].values[0]

                # Determine asterisk type
                if p_val < 0.05 and p_val_log < 0.05:
                        x1, x2 = i * 2, i * 2 + 1
                        y1, y2 = means[x1] + stds[x1], means[x2] + stds[x2]
                        max_y = max(y1, y2)
                        h = 0.05 * max(means)

                        # Check odds ratio confidence interval
                        if ci_lower > 1 or ci_upper < 1:
                                if odds > 1.1 or odds < 0.99 :
                                        star = '**'
                                else:
                                        star = '*'
                        else :
                                star = '*'

                        ax.plot([x1, x1, x2, x2], [max_y+h, max_y+2*h, max_y+2*h, max_y+h], color='black')
                        ax.text((x1+x2)/2, max_y + 2.5*h, star, ha='center', va='bottom', fontsize=14)

        plot_file = os.path.join(save_path, "von_plot.png")
        print(plot_file)
        plt.savefig(plot_file)
        plt.close()

        metrics = ["mean_eval_max", "max_eval_max", "mean_UQ_eval_max"]
        label = ["mean_eval_max_homo", "mean_eval_max_hete", "max_eval_max_homo", "max_eval_max_hete", "mean_UQ_eval_max_homo","mean_UQ_eval_max_hete"]
        means = [dfs[region][field].mean() for field in label]
        stds = [dfs[region][field].std() for field in label]

        fig, ax = plt.subplots(figsize=(7, 7))
        bar_colors = ['#ff9999', '#66b3ff'] * len(metrics)
        ax.bar(label, means, yerr=stds, capsize=7, color=bar_colors)
        ax.set_xticks(range(len(label)))
        ax.set_xticklabels([l.replace("_homo", " (Homo)").replace("_hete", " (Hete)") for l in label], rotation=20)
        ax.set_title(f"Comparison The Maximum Principal Stress in {region}")
        ax.set_ylabel("Von Mises Stress (dyn/cm²)")

        legend_elements = [
                Patch(facecolor='#ff9999', label='Homogeneous'),
                Patch(facecolor='#66b3ff', label='Heterogeneous')
        ]
        ax.legend(handles=legend_elements)

        # Add asterisks if significant
        for i, metric in enumerate(metrics):
                row = stat[(stat["region"] == region) & (stat["field"] == metric)]
                if not row.empty:
                        p_val = row["p_paired_test"].values[0]
                        p_val_log = row["logit_p"].values[0]
                        odds = row["odds_ratio"].values[0]
                        ci_lower = row["odds_ratio_CI_lower"].values[0]
                        ci_upper = row["odds_ratio_CI_upper"].values[0]

                # Determine asterisk type
                if p_val < 0.05 and p_val_log < 0.05:
                        x1, x2 = i * 2, i * 2 + 1
                        y1, y2 = means[x1] + stds[x1], means[x2] + stds[x2]
                        max_y = max(y1, y2)
                        h = 0.05 * max(means)

                        # Check odds ratio confidence interval
                        if ci_lower > 1 or ci_upper < 1:
                                if odds > 1.1 or odds < 0.99 :
                                        star = '**'
                                else:
                                        star = '*'
                        else :
                                star = '*'

                        ax.plot([x1, x1, x2, x2], [max_y+h, max_y+2*h, max_y+2*h, max_y+h], color='black')
                        ax.text((x1+x2)/2, max_y + 2.5*h, star, ha='center', va='bottom', fontsize=14)

        plot_file = os.path.join(save_path, "eval_max.png")
        print(plot_file)
        plt.savefig(plot_file)
        plt.close()

        metrics = ["concen_von", "concen_eval_max"]
        label = ["concen_von_homo", "concen_von_hete", "concen_eval_max_homo", "concen_eval_max_hete"]
        means = [dfs[region][field].mean() for field in label]
        stds = [dfs[region][field].std() for field in label]

        fig, ax = plt.subplots(figsize=(7, 7))
        bar_colors = ['#ff9999', '#66b3ff'] * len(metrics)
        ax.bar(label, means, yerr=stds, capsize=7, color=bar_colors)
        ax.set_xticks(range(len(label)))
        ax.set_xticklabels([l.replace("_homo", " (Homo)").replace("_hete", " (Hete)") for l in label], rotation=20)
        ax.set_title(f"Comparison The Concentration of Stress in {region}")
        ax.set_ylabel("Concentration of Stress")

        legend_elements = [
                Patch(facecolor='#ff9999', label='Homogeneous'),
                Patch(facecolor='#66b3ff', label='Heterogeneous')
        ]
        ax.legend(handles=legend_elements)

        # Add asterisks if significant
        for i, metric in enumerate(metrics):
                row = stat[(stat["region"] == region) & (stat["field"] == metric)]
                if not row.empty:
                        p_val = row["p_paired_test"].values[0]
                        p_val_log = row["logit_p"].values[0]
                        odds = row["odds_ratio"].values[0]
                        ci_lower = row["odds_ratio_CI_lower"].values[0]
                        ci_upper = row["odds_ratio_CI_upper"].values[0]

                # Determine asterisk type
                if p_val < 0.05 and p_val_log < 0.05:
                        x1, x2 = i * 2, i * 2 + 1
                        y1, y2 = means[x1] + stds[x1], means[x2] + stds[x2]
                        max_y = max(y1, y2)
                        h = 0.05 * max(means)

                        # Check odds ratio confidence interval
                        if ci_lower > 1 or ci_upper < 1:
                                if odds > 1.1 or odds < 0.99 :
                                        star = '**'
                                else:
                                        star = '*'
                        else :
                                star = '*'

                        ax.plot([x1, x1, x2, x2], [max_y+h, max_y+2*h, max_y+2*h, max_y+h], color='black')
                        ax.text((x1+x2)/2, max_y + 2.5*h, star, ha='center', va='bottom', fontsize=14)

        plot_file = os.path.join(save_path, "concen_stress.png")
        print(plot_file)
        plt.savefig(plot_file)
        plt.close() 
 
        metrics = ["mean_eval_ratio", "concen_udir"]
        label = ["mean_eval_ratio_homo", "mean_eval_ratio_hete", "concen_udir_homo", "concen_udir_hete"]
        means = [dfs[region][field].mean() for field in label]
        stds = [dfs[region][field].std() for field in label]

        fig, ax = plt.subplots(figsize=(7, 7))
        bar_colors = ['#ff9999', '#66b3ff'] * len(metrics)
        ax.bar(label, means, yerr=stds, capsize=7, color=bar_colors)
        ax.set_xticks(range(len(label)))
        ax.set_xticklabels([l.replace("_homo", " (Homo)").replace("_hete", " (Hete)") for l in label], rotation=20)
        ax.set_title(f"Comparison Directionality of Stress in {region}")
        ax.set_ylabel("Ratio")

        legend_elements = [
                Patch(facecolor='#ff9999', label='Homogeneous'),
                Patch(facecolor='#66b3ff', label='Heterogeneous')
        ]
        ax.legend(handles=legend_elements)

        # Add asterisks if significant
        for i, metric in enumerate(metrics):
                row = stat[(stat["region"] == region) & (stat["field"] == metric)]
                if not row.empty:
                        p_val = row["p_paired_test"].values[0]
                        p_val_log = row["logit_p"].values[0]
                        odds = row["odds_ratio"].values[0]
                        ci_lower = row["odds_ratio_CI_lower"].values[0]
                        ci_upper = row["odds_ratio_CI_upper"].values[0]

                # Determine asterisk type
                if p_val < 0.05 and p_val_log < 0.05:
                        x1, x2 = i * 2, i * 2 + 1
                        y1, y2 = means[x1] + stds[x1], means[x2] + stds[x2]
                        max_y = max(y1, y2)
                        h = 0.05 * max(means)

                        # Check odds ratio confidence interval
                        if ci_lower > 1 or ci_upper < 1:
                                if odds > 1.1 or odds < 0.99 :
                                        star = '**'
                                else:
                                        star = '*'
                        else :
                                star = '*'

                        ax.plot([x1, x1, x2, x2], [max_y+h, max_y+2*h, max_y+2*h, max_y+h], color='black')
                        ax.text((x1+x2)/2, max_y + 2.5*h, star, ha='center', va='bottom', fontsize=14)

        plot_file = os.path.join(save_path, "direc_stress.png")
        print(plot_file)
        plt.savefig(plot_file)
        plt.close() 
