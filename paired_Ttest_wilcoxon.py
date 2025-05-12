import pandas as pd
import numpy as np
import scipy.stats as stats
import pingouin as pg
import statsmodels.api as sm
from tqdm import tqdm
import os
import shutil
import sys
from utils.mystat import *

# Define folders
dir_runfebio = "/dagon1/achitsaz/runfebio/" 
table_dir = "./table/"
save_dir = "./stat/"

list_of_region = ["aneu", "dome", "body", "neck", "thin", "thick", "rupt", "bleb"]
list_of_field = [
    "mean_von", "mean_UQ_von", "max_von", "concen_von",
    "mean_eval_max", "mean_UQ_eval_max", "max_eval_max", "concen_eval_max",
    "force_ratio", "area_ratio"
]
list_of_region_class = ["blind_class", "thick", "thin", "hete"]

list_of_region = ["aneu"]
list_of_field = [
    "mean_von", "mean_UQ_von", "max_von", "concen_von",
    "mean_eval_max", "mean_UQ_eval_max", "max_eval_max", "concen_eval_max",
    "force_ratio", "area_ratio",
    "mean_eval_ratio",
    "concen_udir","force_udir_ratio","area_udir_ratio"
]
list_of_region_class = ["blind_class"]
# Clean save directory
if os.path.exists(save_dir):
    shutil.rmtree(save_dir)
os.makedirs(save_dir, exist_ok=True)

# Load all dataframes
dfs = {region: pd.read_csv(f"{table_dir}{region}/continous_fields.csv") for region in list_of_region}

# Prepare results list
summary_dataset = []

# Loop over all regions and fields
for statusID in range(1):
    if statusID == 0:
        status_name = "blind_status"
    if statusID == 1:
        status_name = "ruptured"  
    if statusID == 2:
        status_name = "stable"      
    for region_class in list_of_region_class:
        results = []
        output_csv = os.path.join(save_dir, f"{status_name}_{region_class}_results.csv")

        for region in list_of_region:
            for field in list_of_field:
                # print(f"Processing: Status={status_name}Region={region}, Field={field}")
                homo_field = []
                hete_field = []
                diff_field = []
                max_odds_ratio = 0
                try:
                    if statusID == 0:    
                        if region_class != "blind_class":
                            filtered = dfs[region][dfs[region]["region_class"] == region_class]
                        else :
                            filtered = dfs[region]
                    if statusID == 1:    
                        if region_class != "blind_class":
                            filtered = dfs[region][
                                (dfs[region]["region_class"] == region_class) & (dfs[region]["status"] == True)
                            ]
                        else :
                            filtered = dfs[region][dfs[region]["status"] == True]     
                    if statusID == 2:    
                        if region_class != "blind_class":
                            filtered = dfs[region][
                                (dfs[region]["region_class"] == region_class) & (dfs[region]["status"] == False)
                            ]
                        else :
                            filtered = dfs[region][dfs[region]["status"] == False]          
                    
                    homo_field = np.array(filtered[f"{field}_homo"])
                    hete_field = np.array(filtered[f"{field}_hete"])
                    diff_field = hete_field - homo_field
                    # check number of observation in the whole table 
                    if region == "aneu":
                        num_obser = len(homo_field)
                    # check possiblity of status*region_class
                    if len(homo_field) == 0:
                        # print(f"\u2718 There is not data which aneurysm is {status_name} within {region_class} group!\n")
                        continue
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
                    
                    # Fitting Logistic Regression
                    # Combine & standardize
                    all_field = np.concatenate([homo_field, hete_field])
                    condition = np.array([0]*len(homo_field) + [1]*len(hete_field))
                    if len(np.unique(condition)) < 2:
                        print(f"Skipping {field}: only one class present", file=sys.stderr)
                        continue

                    # z-score
                    mean_ = all_field.mean()
                    std_  = all_field.std(ddof=0)
                    if std_ == 0:
                        std_ = 1.0
                    all_field_std = (all_field - mean_) / std_

                    # Build logistic dataset
                    df_logit = pd.DataFrame({
                        'stress_value': all_field_std,
                        'condition': condition
                    })

                    # Fit logistic regression with a try/except
                    X = sm.add_constant(df_logit[['stress_value']])
                    y = df_logit['condition']
                    logit_model = sm.Logit(y, X)

                    try:
                        logit_result = logit_model.fit(disp=0)
                    except Exception as e:
                        print(f"Logit failed for {field}: {e}", file=sys.stderr)
                        continue

                    coef = logit_result.params['stress_value']
                    p_logit = logit_result.pvalues['stress_value']
                    # Odds ratio and 95% CI
                    odds_ratio = np.exp(coef)
                    ci_lower, ci_upper = np.exp(logit_result.conf_int().loc['stress_value'])
                    # Bootstrapt Logistic Regressin 
                    n_bootstrap = 1000
                    boot_odds_ratios = []

                    for _ in tqdm(range(n_bootstrap)):
                        # Resample indices with replacement
                        sample_idx = np.random.choice(len(df_logit), len(df_logit), replace=True)
                        sample_df = df_logit.iloc[sample_idx]
                        
                        X_boot = sm.add_constant(sample_df[['stress_value']])
                        y_boot = sample_df['condition']
                        
                        try:
                            model = sm.Logit(y_boot, X_boot).fit(disp=0)
                            coef_boot = model.params['stress_value']
                            or_boot = np.exp(coef_boot)
                            boot_odds_ratios.append(or_boot)
                        except:
                            continue  # skip iterations where model fails to converge

                    # Compute bootstrap statistics
                    boot_odds_ratios = np.array(boot_odds_ratios)
                    median_or = np.median(boot_odds_ratios)
                    boot_ci_lower = np.percentile(boot_odds_ratios, 2.5)
                    boot_ci_upper = np.percentile(boot_odds_ratios, 97.5)
                    
                    # Append result
                    results.append({
                        "region": region,
                        "field": field,
                        "num_obser" : len(homo_field),
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
                        "odds_ratio_CI_upper": ci_upper,
                        "n_bootstrap": n_bootstrap,
                        "boot_median_or": median_or,
                        "boot_ci_lower": boot_ci_lower,
                        "boot_ci_upper": boot_ci_upper

                    })
                except Exception as e:
                    # print(f"Failed for Status={status_name},class={region_class}, region={region}, field={field} due to {e}", file=sys.stderr)
                    continue 
        # check possiblity of status*region_class
        if num_obser == 0:
            # print(f"\u2718 There is not data which aneurysm is {status_name} within {region_class} group!\n")
            continue
        # Save all results of statistics into one CSV
        df_results = pd.DataFrame(results)
        df_results.to_csv(output_csv, index=False)
        # print(f"\n✅ All results saved to: {output_csv}")

        # Define thresholds
        paired_threshold = 0.05
        logit_threshold = 0.05
        # Find significant results
        significant = df_results[(df_results["p_paired_test"] < paired_threshold) | (df_results["logit_p"] < logit_threshold)]
        # significant = df_results[(df_results["p_paired_test"] < paired_threshold) ]
        significant = df_results
        # Sort them nicely
        significant = significant.sort_values(by=["region", "field"])
        # Show details
        print(f"\n✅-----------------------Show Details of {region_class} Cases --------------------------------- \n")
        for idx, row in significant.iterrows():
            region = row["region"]
            field = row["field"]
            num_obser = row["num_obser"]
            
            mean_homo = row["mean_homo"]
            mean_hete = row["mean_hete"]
            p_paired = row["p_paired_test"]
            stat_used = row["paired_test_used"]
            
            logit_p = row["logit_p"]
            odds_ratio = row["odds_ratio"]
            ci_lower = row["odds_ratio_CI_lower"]
            ci_upper = row["odds_ratio_CI_upper"]

            boot_odds_ratio = row["boot_median_or"]
            boot_ci_lower = row["boot_ci_lower"]
            boot_ci_upper = row["boot_ci_upper"]
            n_bootstrap = row["n_bootstrap"]


            # if region != "aneu":
            #     continue

            # if ci_lower < 1 and ci_upper > 1 :
            #     continue
            print(f"---")
            print(f"📌 status: {status_name} | class: {region_class} | Region: {region} | Field: {field} | num_obser: {num_obser}")
            print(f"    - Paired test used: {stat_used}")
            print(f"    - Mean Homo: {mean_homo:.4f} | Mean Hete: {mean_hete:.4f}")
            print(f"    - p-value (Paired Test): {p_paired:.4e}")
            print(f"    - p-value (Logistic Regression): {logit_p:.4e}")
            print(f"    - Odds Ratio: {odds_ratio:.4f} (95% CI: [{ci_lower:.4f}, {ci_upper:.4f}])")
            print(f"    - Median Odds Ratios (Bootstrapt {n_bootstrap}) {boot_odds_ratio:.4f} (95% CI: [{boot_ci_lower:.4f}, {boot_ci_upper:.4f}])")
            
            # Interpretation
            if odds_ratio > 1:
                interpretation = "higher odds for heterogeneous condition"
            else:
                interpretation = "lower odds for heterogeneous condition"
            print(f"    - Interpretation: {interpretation}")
            print("")



