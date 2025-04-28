import numpy as np

def ordinary_stat (data):
    mean = np.mean(data)
    std = np.std(data, ddof=0)
    median = np.median(data)
    q1 = np.percentile(data, 25)  # 25th percentile (Q1)
    q3 = np.percentile(data, 75)  # 75th percentile (Q3)
    iqr = q3 - q1
    return mean, std, median, iqr 