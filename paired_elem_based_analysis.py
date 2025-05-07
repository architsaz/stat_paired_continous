import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
import os
import shutil
import sys

dir_runfebio = "/dagon1/achitsaz/runfebio/" 
save_dir = "./stat.elem/"
studies = ["pst.1", "pst.2"]
region = "aneu"
all_results = []
dir_approved_case = dir_runfebio+"successful_cases.txt"
with open(dir_approved_case,"r") as f:
    list_approved_case = [line.strip() for line in f]
# clean save data directory    
shutil.rmtree(save_dir)
os.makedirs(save_dir, exist_ok=True)      
for case in list_approved_case:   
    region_table = f'{save_dir}{case}'
    if os.path.exists(region_table):
        shutil.rmtree(region_table)
    os.makedirs(region_table, exist_ok=True)      
# list_approved_case = ["a06161.1","a10151.1", "agh101.4"]
for case in list_approved_case:
    # global arraies for each case 
    accepted_case = True
    npoin = 0
    nelem = 0
    region_p = []
    elems = []
    pts = []
    color_mask = []
    bleb_mask = []
    homo_von = []
    hete_von = []
    homo_eval_max = []
    hete_eval_max = []
    homo_eval_min = []
    hete_eval_min = []    
    homo_eval_ratio = []
    hete_eval_ratio = []   
    area_region = [] 

    # Section 1 : read all required fields 
    try:
        find_region = False
        find_elems = False
        find_point = False
        path_file = f'{dir_runfebio}{case}/data/labels_srf.zfem'
        lines_to_skip = 0
        with open (path_file,"r") as f:
            for line in f:
                line = line.strip()
                if find_point == True and lines_to_skip == 1 :
                    npoin = int(line)
                if find_elems == True and lines_to_skip == 1 :
                    nelem = int(line)    
                if lines_to_skip > 0:
                    lines_to_skip -= 1
                    continue
                if "regions" in line:
                    find_region = True
                    lines_to_skip = 1
                    continue
                if find_region:
                    if "ELEMENTS" in line:
                        find_region = False
                        lines_to_skip = 0
                        continue
                if "TRIANGLE" in line:
                    find_elems = True
                    lines_to_skip = 2
                    continue
                if find_elems:
                    if "END" in line:
                        find_elems = False
                        lines_to_skip = 0
                        continue    
                if "POINTS" in line:
                    find_point = True
                    lines_to_skip = 1
                    continue
                if find_point:
                    if "REAL" in line:
                        find_point = False
                        lines_to_skip = 0
                        continue
                if lines_to_skip == 0 :
                    if find_point:
                        value = line.split()
                        if len(value) == 3:
                            pts.append([float(v) for v in value])
                        else:
                            print(f'ERROR: number of element in {case} line : {line} is not 3', file=sys.stderr)
                            accepted_case = False
                            continue

                    if find_elems:
                        value = line.split()
                        if len(value) == 3:
                            elems.append([int(v) for v in value])
                        else:
                            print (f'ERROR: number of element in {case} line : {line} is not 3', file=sys.stderr)
                            accepted_case = False
                            continue
                        
                    if find_region:
                        value = line.split()
                        if len(value) == 1:
                            region_p.append(int(value[0]))
                        else:
                            print (f'ERROR: number of element in {case} line : {line} is not 1', file=sys.stderr)
                            accepted_case = False
                            continue

            if len(region_p) != npoin:
                print (f'ERROR: can not read regions field in case {case}!', file =sys.stderr)
                accepted_case = False
                continue

            if len(pts) != npoin:
                print (f'ERROR: can not read pts field in case {case}!', file =sys.stderr)
                accepted_case = False
                continue

            if len(elems) != nelem:
                print (f'ERROR: can not read ELEMENTS field in case {case}!', file =sys.stderr)
                accepted_case = False
                continue
    except FileExistsError:
        print(f"ERROR: labels_surf.zfem does not exist in the data directory of case {case}", file=sys.stderr)
        accepted_case = False
        continue
    # check the status of case
    if accepted_case == False :
        continue      
    # convert region_p to region_e
    region_e = []
    for elem in elems: 
        pointID = elem[0]-1
        region_e.append(region_p[pointID])    
    # save color mask in color_mask
    try:
        file_path=dir_runfebio+case+f"/data/{case}.wall"
        with open(file_path,"r") as file:
            for line in file:
                row = line.strip()
                value = row.split()
                if (len(value)==1):
                    color_mask.append(int(value[0]))
                else:
                    print(f"ERROR in reading color mask: number of element in line : {row} is more than 1. case{case}", file=sys.stderr)
                    accepted_case = False
                    continue 
    except FileExistsError:
        print(f"ERROR: color mask does not exist in the data directory of case {case}", file=sys.stderr)
        accepted_case = False
        continue
    if len(color_mask) != nelem:
        print(f"ERROR: number of element in color mask is not match with nelem", file=sys.stderr)
        accepted_case = False
        continue
    # check the status of case
    if accepted_case == False :
        continue
    # check existance at least one thick or thin color region
    exist_one_color = False
    for ele in range(len(color_mask)):
        if color_mask[ele] in [1,4,7]:
            exist_one_color = True
    if exist_one_color == False:
        print(f"ERROR in color mask: there is no thin or thick region in the color mask of case {case}", file=sys.stderr) 
        continue
    # find bleb mask 
    path_file = f'{dir_runfebio}{case}/data/{case}.bleb'
    if os.path.exists(path_file):
        try:
            with open (path_file,"r") as file:
                for line in file:
                    row = line.strip()
                    value = row.split()
                    if (len(value)==1):
                        bleb_mask.append(int(value[0]))
                    else:
                        print(f"ERROR in reading bleb mask: number of element in line : {row} is more than 1. case{case}", file=sys.stderr)
                        accepted_case = False
                        continue 
        except Exception as e:
            print(f"An error occurred: {e}")
        # check number of element in bleb_mask field
        if len (bleb_mask) != nelem :
            print(f"ERROR: number of element in .beleb mask for {case} is not match with nelem!", file=sys.stderr)
            accepted_case = False
            continue     
    else:
        print(f"The file '{path_file}' does not exist for case {case}", file=sys.stderr)    
    # Section 2 : read all categorial data fields for each condition   
    for study in studies:
        try:
            path_file = f'{dir_runfebio}{case}/{study}/stress_analysis_0.vtk'
            von_mises = []
            eval_max = []
            eval_min = []
            eval_ratio = []
            area = []
            mode = None
            lines_to_skip = 0
            num_fields = 0
            with open(path_file, "r") as f:
                for line in f:
                    line = line.strip()

                    if lines_to_skip > 0:
                        lines_to_skip -= 1
                        continue

                    if "area" in line:
                        mode = "area"
                        lines_to_skip = 2
                        num_fields += 1
                        continue

                    if "Shear_Von_Mises_Stress" in line:
                        mode = "von_mises"
                        lines_to_skip = 2
                        num_fields += 1
                        continue

                    if "EValue_max" in line:
                        mode = "eval_max"
                        lines_to_skip = 2
                        num_fields += 1
                        continue

                    if "EValue_min" in line:
                        mode = "eval_min"
                        lines_to_skip = 2
                        num_fields += 1
                        continue

                    if "Eval_ratio" in line:
                        mode = "eval_ratio"
                        lines_to_skip = 2
                        num_fields += 1
                        continue

                    if "SCALARS" in line:
                        mode = None
                        lines_to_skip = 0
                        continue

                    values = line.split()
                    if mode != None:
                        if len(values) == 1: 
                            if mode == "area":
                                area.append(float(values[0]))  # or use the right index
                            elif mode == "von_mises":
                                von_mises.append(abs(float(values[0]))) 
                            elif mode == "eval_max":
                                eval_max.append(abs(float(values[0]))) 
                            elif mode == "eval_min":
                                eval_min.append(abs(float(values[0])))                                  
                            elif mode == "eval_ratio":      
                                eval_ratio.append(float(values[0]))
                        else:
                            print(f'ERROR: number of element in {line} is more than expected!', file=sys.stderr)
                            accepted_case = False
                            continue        
                if num_fields !=  5:
                    print(f'ERROR: can not find all field in case {case} study {study}!', file=sys.stderr)
                    accepted_case = False
                    continue 
                if len(area) != nelem:
                    print (f'ERROR: can not read Area field in case {case} study {study}!', file=sys.stderr)
                    accepted_case = False
                    continue 
                if len(von_mises) != nelem:
                    print (f'ERROR: can not read von_mises field in case {case} study {study}!', file=sys.stderr)
                    accepted_case = False
                    continue    
                if len(eval_max) != nelem:
                    print (f'ERROR: can not read eval_max field in case {case} study {study}!', file=sys.stderr)
                    accepted_case = False
                    continue
                if len(eval_ratio) != nelem:
                    print (f'ERROR: can not read eval_ratio field in case {case} study {study}!', file=sys.stderr)
                    accepted_case = False
                    continue             
        except FileExistsError:
            print(f"ERROR: stress_analysis does not exist in the pst.# directory of case {case}", file=sys.stderr)
            accepted_case = False
            continue            
        #check status of case 
        if accepted_case == False :
            continue  
        # make full_stress_cls
        if study == "pst.1":
            for i in range(nelem):
                if region_e[i] in [16,8,4]:    
                    homo_von.append(von_mises[i])
                    homo_eval_max.append(eval_max[i])
                    homo_eval_min.append(eval_min[i])
                    homo_eval_ratio.append(eval_ratio[i])
                    area_region.append(area[i])
        if study == "pst.2":
            for i in range(nelem):
                if region_e[i] in [16,8,4]:    
                    hete_von.append(von_mises[i])
                    hete_eval_max.append(eval_max[i])
                    hete_eval_min.append(eval_min[i])
                    hete_eval_ratio.append(eval_ratio[i])        
    #check status of case 
    if accepted_case == False :
        continue 
    for i in range(8):
        if i == 0:
            homo_field = homo_von
            hete_field = hete_von
            field_name = "Von Mises Stress"
            field_name_abb = "von_misses"
        if i == 1:
            homo_field = homo_eval_max
            hete_field = hete_eval_max
            field_name = "Max of Principal Stress"
            field_name_abb = "eval_max"
        if i == 2:
            homo_field = homo_eval_min
            hete_field = hete_eval_min
            field_name = "Min of Principal Stress"
            field_name_abb = "eval_min"    
        if i == 3:
            homo_field = homo_eval_ratio
            hete_field = hete_eval_ratio   
            field_name = "Eigenvalue Ratio"
            field_name_abb = "eval_ratio"
        if i == 4:
            homo_field = np.array(homo_von)*np.array(area_region)
            hete_field = np.array(hete_von)*np.array(area_region)
            field_name = "Force based on Von Mises"
            field_name_abb = "fvon_misses"
        if i == 5:
            homo_field = np.array(homo_eval_max)*np.array(area_region)
            hete_field = np.array(hete_eval_max)*np.array(area_region)
            field_name = "force based on Max of Principal Stress"
            field_name_abb = "feval_max"
        if i == 6:
            homo_field = np.array(homo_eval_min)*np.array(area_region)
            hete_field = np.array(hete_eval_min)*np.array(area_region)
            field_name = "force based on Min of Principal Stress"
            field_name_abb = "feval_min"            
        if i == 7:
            homo_field = np.array(homo_eval_ratio)*np.array(area_region)
            hete_field = np.array(hete_eval_ratio)*np.array(area_region)  
            field_name = "area weighted Eigenvalue Ratio"
            field_name_abb = "AWeval_ratio"                  
        # Create DataFrame
        df = pd.DataFrame({
            'homogeneous': homo_field,
            'heterogeneous': hete_field
        })

        # Fit linear regression model
        X = df[['homogeneous']]  # Independent variable
        y = df['heterogeneous']  # Dependent variable
        model = LinearRegression().fit(X, y)
        # Get slope, intercept, R²
        slope = model.coef_[0]
        intercept = model.intercept_
        r2 = model.score(X, y)

        # Save regression results to CSV
        results_df = pd.DataFrame({
            "case": case,
            "field": field_name,
            "region": region,
            "slope": [slope],
            "intercept": [intercept],
            "r_squared": [r2]
        })
        all_results.append(results_df)

        # Create x values for the regression line
        x_vals = pd.DataFrame(np.linspace(X.min().iloc[0], X.max().iloc[0], 100), columns=["homogeneous"])
        y_vals = model.predict(x_vals)

        # Plot
        plt.figure(figsize=(8, 6))
        plt.scatter(X, y, label='Data points')
        plt.plot(x_vals, y_vals, color='red', label='Regression line')
        plt.xlabel(f"Homogeneous {field_name}")
        plt.ylabel(f"Heterogeneous {field_name}")
        plt.title("Linear Regression: Homogeneous vs Heterogeneous")

        # Add annotation
        textstr = f"Slope: {slope:.3f}\nIntercept: {intercept:.3f}\nR²: {r2:.3f}"
        plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
                fontsize=10, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))

        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{save_dir}{case}/{field_name_abb}.png', dpi=300)
        plt.close()  

# Concatenate all individual DataFrames into one
final_results = pd.concat(all_results, ignore_index=True)
# Group by 'field' and compute mean and std
summary_stats = final_results.groupby("field").agg({
    "slope": ['mean', 'std'],
    "intercept": ['mean', 'std'],
    "r_squared": ['mean', 'std']
}).reset_index()
# Optional: flatten column names
summary_stats.columns = ['field', 'slope_mean', 'slope_std',
                         'intercept_mean', 'intercept_std',
                         'r_squared_mean', 'r_squared_std']

# Save summary to CSV
summary_stats.to_csv(f'{save_dir}/LR_summary_stats.csv', index=False)
# Optional: print it
print(summary_stats)
# Save to CSV once after the loop
final_results.to_csv(f'{save_dir}/LR_results.csv', index=False)