
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def sizes():
    """ Create sizes data achieved from the FADN public database set for the year of 2016 in Walloon region

    :return:                number of farmers in each size category
    :type:                  Numpy Array

    """

    # size_cat = ["25000 to 50000", "50000 to 100000", "100000 to 500000", "more than 500000"]
    size_cat = ["XSmall", "Small", "Medium", "Large"]
    size_num = [240, 2960, 6510, 860]
    # size_ha = np.array([1, 16.8, 41.8, 95.6]) # in FADN xsmall value is not give, here assume to be 1 ha

    # size_r_meter = np.sqrt(size_ha * 10000) # 边长
    # size_r_average = np.sqrt(382900) # average ha of arable area given by FADN

    sizes = pd.DataFrame(0, index=np.arange(len(size_cat)), columns=('category', 'number'))
    sizes['category'] = size_cat.copy()
    sizes['number'] = size_num.copy()

    return sizes


def make_size_cdf(var_array):
    
    ser_full = np.zeros(0)

    ser0 = np.zeros(var_array['number'][0])
    ser1 = np.ones(var_array['number'][1])
    ser2 = np.ones(var_array['number'][2]) + 1
    ser3 = np.ones(var_array['number'][3]) + 2
    ser_full = np.append(ser_full, ser0)
    ser_full = np.append(ser_full, ser1)
    ser_full = np.append(ser_full, ser2)
    ser_full = np.append(ser_full, ser3)
        
    h, x1 = np.histogram(ser_full, bins=4, density=True)
    dx = x1[2] - x1[1]
    f1 = np.cumsum(h) * dx
    perc = np.column_stack(([0, 1, 2, 3], f1))
    
    return perc
    

    
    

def farmer_data(size_cdf, switch_size):

    ts = np.random.random_sample()
                
    tt = np.where(size_cdf[:, [1]] >= ts)
    ten_stat = min(tt[0])


    if ten_stat == size_cdf[0, [0]]:  # XS farm, switching averse. k_s: k for size
        k_s = 0  
    elif ten_stat == size_cdf[1, [0]]:  # S farm, switching neutral
        k_s = 2
    elif ten_stat == size_cdf[2, [0]]:  # M farm, switching neutral
        k_s = 2
    elif ten_stat == size_cdf[3, [0]]:  # L, switching tolerant
        k_s = 1

    
    alpha_s = switch_size[k_s][0]
    beta_s = switch_size[k_s][1]

    agent_data = {
            "FarmSize": ten_stat,
            "SizeAlpha": alpha_s,
            "SizeBeta": beta_s,
            "nFields": 1
                }
    return agent_data



def company_data(lc):


    if lc == 6:
        product = 'biomethane'
    if lc == 7:
        product = 'combined heat and power'
    if lc == 8:
        product = 'bioetanol'

    agent_data = {"Product": product}

    return agent_data

    