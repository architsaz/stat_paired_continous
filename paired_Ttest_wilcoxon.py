import pandas as pd
import numpy as np
import scipy.stats as stats
import pingouin as pg
import statsmodels.api as sm
import os
import shutil
import sys
from utils.mystat import *

# Define folders
dir_runfebio = "/dagon1/achitsaz/runfebio/" 
table_dir = "./table/"
save_dir = "./stat/"
output_csv = os.path.join(save_dir, "summary_results.csv")

list_of_region = ["aneu", "dome", "body", "neck", "thin", "thick", "rupt", "bleb"]
list_of_field = [
    "mean_von", "mean_UQ_von", "max_von", "concen_von",
    "mean_eval_max", "mean_UQ_eval_max", "max_eval_max", "concen_eval_max"
]

# Clean save directory
if os.path.exists(save_dir):
    shutil.rmtree(save_dir)
os.makedirs(save_dir, exist_ok=True)

# Load all dataframes
dfs = {region: pd.read_csv(f"{table_dir}{region}/continous_fields.csv") for region in list_of_region}

# Prepare results list
results = []

# Loop over all regions and fields
for region in list_of_region:
    for field in list_of_field:
        print(f"Processing: Region={region}, Field={field}")
        
        try:
            homo_field = np.array(dfs[region][f"{field}_homo"])
            hete_field = np.array(dfs[region][f"{field}_hete"])
            diff_field = hete_field - homo_field

            # Sanity check
            if len(homo_field) != len(hete_field):
                print(f"ERROR: Size mismatch for {region}-{field}!", file=sys.stderr)
                continue
            if np.allclose(diff_field, 0):
                print(f"WARNING: All differences are zero for {region}-{field}.", file=sys.stderr)
            
            # Descriptive stats
            mean_homo, std_homo, median_homo, iqr_homo = ordinary_stat(homo_field)
            mean_hete, std_hete, median_hete, iqr_hete = ordinary_stat(hete_field)
            
            # Normality
            _, shapiro_p_homo = stats.shapiro(homo_field)
            _, shapiro_p_hete = stats.shapiro(hete_field)
            _, p_diff = stats.shapiro(diff_field)

            # Paired test
            if p_diff > 0.05 and shapiro_p_homo > 0.05 and shapiro_p_hete > 0.05 :
                stat_test = "paired t-test"
                test_result = pg.ttest(homo_field, hete_field, paired=True)
                p_paired = test_result["p-val"].values[0]
            else:
                stat_test = "Wilcoxon signed-rank"
                test_result = pg.wilcoxon(homo_field, hete_field)
                p_paired = test_result["p-val"].values[0]
            
            # Logistic regression
            all_field = np.concatenate([homo_field, hete_field])
            condition = np.array([0]*len(homo_field) + [1]*len(hete_field))
            
            df_logit = pd.DataFrame({
                'stress_value': all_field,
                'condition': condition
            })
            X = sm.add_constant(df_logit[['stress_value']])
            y = df_logit['condition']
            logit_model = sm.Logit(y, X)
            logit_result = logit_model.fit(disp=0)  # silent fit

            coef = logit_result.params['stress_value']
            p_logit = logit_result.pvalues['stress_value']
            # Odds ratio and 95% CI
            odds_ratio = np.exp(coef)
            ci_lower, ci_upper = np.exp(logit_result.conf_int().loc['stress_value'])
            
            # Append result
            results.append({
                "region": region,
                "field": field,
                "mean_homo": mean_homo,
                "std_homo": std_homo,
                "median_homo": median_homo,
                "iqr_homo": iqr_homo,
                "mean_hete": mean_hete,
                "std_hete": std_hete,
                "median_hete": median_hete,
                "iqr_hete": iqr_hete,
                "shapiro_p_homo": shapiro_p_homo,
                "shapiro_p_hete": shapiro_p_hete,
                "shapiro_p_diff": p_diff,
                "paired_test_used": stat_test,
                "p_paired_test": p_paired,
                "logit_coef": coef,
                "logit_p": p_logit,
                "odds_ratio": odds_ratio,
                "odds_ratio_CI_lower": ci_lower,
                "odds_ratio_CI_upper": ci_upper
            })

        except Exception as e:
            print(f"Failed for region={region}, field={field} due to {e}", file=sys.stderr)
            continue

# Save all results into one CSV
df_results = pd.DataFrame(results)
df_results.to_csv(output_csv, index=False)

print(f"\n✅ All results saved to: {output_csv}")

# Load the result csv
df = pd.read_csv("./stat/summary_results.csv")

# Define thresholds
paired_threshold = 0.05
logit_threshold = 0.05

# Find significant results
significant = df[(df["p_paired_test"] < paired_threshold) | (df["logit_p"] < logit_threshold)]

# Sort them nicely
significant = significant.sort_values(by=["region", "field"])

# Show details
for idx, row in significant.iterrows():
    region = row["region"]
    field = row["field"]
    
    mean_homo = row["mean_homo"]
    mean_hete = row["mean_hete"]
    p_paired = row["p_paired_test"]
    stat_used = row["paired_test_used"]
    
    logit_p = row["logit_p"]
    odds_ratio = row["odds_ratio"]
    ci_lower = row["odds_ratio_CI_lower"]
    ci_upper = row["odds_ratio_CI_upper"]
    
    print(f"---")
    print(f"📌 Region: {region} | Field: {field}")
    print(f"    - Paired test used: {stat_used}")
    print(f"    - Mean Homo: {mean_homo:.4f} | Mean Hete: {mean_hete:.4f}")
    print(f"    - p-value (Paired Test): {p_paired:.4e}")
    print(f"    - p-value (Logistic Regression): {logit_p:.4e}")
    print(f"    - Odds Ratio: {odds_ratio:.4f} (95% CI: [{ci_lower:.4f}, {ci_upper:.4f}])")
    
    # Interpretation
    if odds_ratio > 1:
        interpretation = "higher odds for heterogeneous condition"
    else:
        interpretation = "lower odds for heterogeneous condition"
    print(f"    - Interpretation: {interpretation}")
    print("")