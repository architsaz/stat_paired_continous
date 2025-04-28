import numpy as np
import pandas as pd
import os
import shutil
import sys

dir_runfebio = "/dagon1/achitsaz/runfebio/" 
table_dir = "./table/"
studies = ["pst.1", "pst.2"]
list_of_region = ["aneu", "dome", "body", "neck", "thin", "thick", "rupt", "bleb"]
for region in list_of_region:
    # clean save data directory   
    region_table = f'{table_dir}{region}'
    if os.path.exists(region_table):
        shutil.rmtree(region_table)
    os.makedirs(region_table, exist_ok=True) 
dir_approved_case = dir_runfebio+"successful_cases.txt"
with open(dir_approved_case,"r") as f:
    list_approved_case = [line.strip() for line in f]  
# list_approved_case = ["a06161.1", "agh101.4"]
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
                                von_mises.append(float(values[0])) 
                            elif mode == "eval_max":
                                eval_max.append(float(values[0]))       
                        else:
                            print(f'ERROR: number of element in {line} is more than expected!', file=sys.stderr)
                            accepted_case = False
                            continue        
                if num_fields != 3 :
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
                homo_von.append(int(von_mises[i]))
                homo_eval_max.append(int(eval_max[i]))
        if study == "pst.2":
            for i in range(nelem):
                hete_von.append(int(von_mises[i]))
                hete_eval_max.append(int(eval_max[i]))
            #check status of case 
    if accepted_case == False :
        continue 
    # Section 3 : Saved Table for Paired Continous 
    # calculate the continous table for each region
    for region in list_of_region:
        if region == "aneu":
            reg = [16, 8, 4] # aneurysm
            used_region_mask = True
        if region == "dome":    
            reg = [16] # dome
            used_region_mask = True
        if region == "body":    
            reg = [8] # body
            used_region_mask = True
        if region == "neck":
            reg = [4] # neck
            used_region_mask = True
        if region == "thin":    
            reg = [1] # red (thinner region)
            used_region_mask = False
        if region == "thick":     
            reg = [4, 7] # yellow & with thicker region
            used_region_mask = False
        if region == "rupt":    
            reg = [9] # rupture
            used_region_mask = False
        # define regional mask
        if used_region_mask:
            mask = region_e
        else:
            mask = color_mask
        # bleb mask
        if region == "bleb":
            if len(bleb_mask) == nelem :
                mask = bleb_mask 
                reg = [1,2,3,4,5,6]
            elif len(bleb_mask) == 0:
                continue        
            else:
                print(f"ERROR: there is bleb mask but can not read it well!", file=sys.stderr)
                continue
        # check existance of region
        region_exist = False
        for ele in range(len(elems)):
            if mask[ele] in reg: 
                region_exist = True
        if region_exist == False:
            print(f"ERROR: there is no {region} region in {case}", file=sys.stderr)
            continue
        # calculate global statistic parameter for each continous fields 
        mean_von_homo = mean_von_hete = mean_eval_max_homo = mean_eval_max_hete = total_area = 0
        mean_UQ_von_homo = mean_UQ_von_hete  = mean_UQ_eval_max_homo = mean_UQ_eval_max_hete = 0
        area_UQ_von_homo = area_UQ_von_hete  = area_UQ_eval_max_homo = area_UQ_eval_max_hete = 0
        homo_von_region = []
        hete_von_region = []
        homo_eval_max_region = []
        hete_eval_max_region = []
        for ele in range (nelem):
            if mask[ele] in reg:
                homo_von_region.append(homo_von[ele])
                hete_von_region.append(hete_von[ele])
                homo_eval_max_region.append(homo_eval_max[ele])
                hete_eval_max_region.append(hete_eval_max[ele])
        Q3_von_homo = np.quantile(homo_von_region, 0.75)
        Q3_von_hete = np.quantile(hete_von_region, 0.75)
        Q3_eval_max_homo = np.quantile(homo_eval_max_region, 0.75)
        Q3_eval_max_hete = np.quantile(hete_eval_max_region, 0.75)
        for ele in range (nelem):
            if mask[ele] in reg: 
                mean_von_homo += homo_von[ele]*area[ele]
                mean_von_hete += hete_von[ele]*area[ele]
                mean_eval_max_homo += homo_eval_max[ele]*area[ele]
                mean_eval_max_hete += hete_eval_max[ele]*area[ele]
                total_area += area[ele]
                if homo_von [ele] > Q3_von_homo:
                    mean_UQ_von_homo += homo_von[ele]*area[ele]
                    area_UQ_von_homo += area [ele]
                if hete_von [ele] > Q3_von_hete:
                    mean_UQ_von_hete += hete_von[ele]*area[ele]
                    area_UQ_von_hete += area [ele]    
                if homo_eval_max [ele] > Q3_eval_max_homo:
                    mean_UQ_eval_max_homo += homo_eval_max[ele]*area[ele]
                    area_UQ_eval_max_homo += area [ele]  
                if hete_eval_max [ele] > Q3_eval_max_hete:
                    mean_UQ_eval_max_hete += hete_eval_max[ele]*area[ele]
                    area_UQ_eval_max_hete += area [ele]   
        area_ratio_homo = (area_UQ_von_homo/total_area)
        area_ratio_hete = (area_UQ_von_hete/total_area)    
        force_ratio_homo = (mean_UQ_eval_max_homo/mean_eval_max_homo)   
        force_ratio_hete = (mean_UQ_eval_max_hete/mean_eval_max_hete)                      
        concen_von_homo = (mean_UQ_von_homo/mean_von_homo) * (total_area/area_UQ_von_homo)
        concen_von_hete = (mean_UQ_von_hete/mean_von_hete) * (total_area/area_UQ_von_hete)
        concen_eval_max_homo = (mean_UQ_eval_max_homo/mean_eval_max_homo) * (total_area/area_UQ_eval_max_homo)
        concen_eval_max_hete = (mean_UQ_eval_max_hete/mean_eval_max_hete) * (total_area/area_UQ_eval_max_hete)
        mean_von_homo = mean_von_homo/total_area
        mean_von_hete = mean_von_hete/total_area
        mean_eval_max_homo = mean_eval_max_homo/total_area
        mean_eval_max_hete = mean_eval_max_hete/total_area
        mean_UQ_von_homo = mean_UQ_von_homo/area_UQ_von_homo
        mean_UQ_von_hete  = mean_UQ_von_hete/area_UQ_von_hete
        mean_UQ_eval_max_homo = mean_UQ_eval_max_homo/area_UQ_eval_max_homo
        mean_UQ_eval_max_hete = mean_UQ_eval_max_hete/area_UQ_eval_max_hete
        
        # Save report CSV
        report_data = {
            "case": [case],
            "mean_von_homo": [mean_von_homo],
            "mean_von_hete": [mean_von_hete],
            "mean_eval_max_homo": [mean_eval_max_homo],
            "mean_eval_max_hete": [mean_eval_max_hete],
            "max_von_homo": [max(homo_von)],
            "max_von_hete": [max(hete_von)],
            "max_eval_max_homo": [max(homo_eval_max)],
            "max_eval_max_hete": [max(hete_eval_max)],
            "mean_UQ_von_homo": [mean_UQ_von_homo],
            "mean_UQ_von_hete": [mean_UQ_von_hete],
            "mean_UQ_eval_max_homo": [mean_UQ_eval_max_homo],
            "mean_UQ_eval_max_hete": [mean_UQ_eval_max_hete],
            "concen_von_homo": [concen_von_homo],
            "concen_von_hete": [concen_von_hete],
            "concen_eval_max_homo": [concen_eval_max_homo],
            "concen_eval_max_hete": [concen_eval_max_hete],
            "force_ratio_homo": [force_ratio_homo],
            "force_ratio_hete": [force_ratio_hete],
            "area_ratio_homo": [area_ratio_homo],
            "area_ratio_hete": [area_ratio_hete]
        }
        region_table = f'{table_dir}{region}'
        report_file = os.path.join(region_table, "continous_fields.csv")
        # Check if the file already exists
        write_header = not os.path.exists(report_file)
        # Append to the file
        pd.DataFrame(report_data).to_csv(report_file, mode='a', index=False, header=write_header)




