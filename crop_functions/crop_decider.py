import numpy as np
import scipy.special as sp
from numpy.lib.stride_tricks import as_strided

import sys
sys.path.append('/Users/tianran/Desktop/dtlca')

import initialize_agents_domain as init_agent
import matplotlib.pyplot as plt


def define_seed(seed):
    """ Creates seed for random selection for testing
    """
    global seed_val

    seed_val = seed

    return


def neighbors(mat, row, col, radius=1):

    rows, cols = len(mat), len(mat[0])
    out = []

    for i in range(row - radius, row + radius + 1):
        row = []
        for j in range(col - radius, col + radius + 1):

            if 0 <= i < rows and 0 <= j < cols:
                row.append(mat[i][j])
            else:
                row.append(0)

        out.append(row)
        
    return np.array(out)


def switching_prob_curve(alpha, beta, n):
    """ Creates probability curves that show likelihood of switching crops based on profits
    """
    x = np.linspace(0, 1.0, num=n)

    fx = sp.betainc(alpha, beta, x)

    return fx


def decide2switch(SizeAlpha, SizeBeta, FamiliarMiscAlpha, FamiliarMiscBeta, fmin, fmax, n, profit, profit_p, profit_neighbor, seed=True):
    if seed:

        try:
            seed_val
        except NameError:
            print("Random seed needs to be initialized using the CropDecider.DefineSeed() Function")

        np.random.seed(seed_val)


    if (profit_neighbor > profit) and (profit_p > profit):

        x = np.linspace(fmin * profit, fmax * profit, num=n)
       
        fx_s = switching_prob_curve(SizeAlpha, SizeBeta, n) # probability of switching with alpha beta for size data
        prob_switch_misc_s = np.interp((profit_p+profit_neighbor)/2, x, fx_s)

        fx_f_misc = switching_prob_curve(FamiliarMiscAlpha, FamiliarMiscBeta, n) # familiarity
        prob_switch_misc_f = np.interp((profit_p+profit_neighbor)/2, x, fx_f_misc)

        prob_switch_misc = prob_switch_misc_s * prob_switch_misc_f

        rand_val = np.random.rand(1)

        if (rand_val < prob_switch_misc):  
            return 1  # Switch
        else:
            return 0  # Do not switch

    else:
        return 0  # Do not switch if not profitable

def assess_profit_potential(crop, profits_current, profit_signals, yield_avg, location_yield, num_crops, crop_ids):
    # Existing Crop ID
    cur_crop_choice_ind = crop.astype('int')

    # assess current and future profit of that given crop
    if np.isin(cur_crop_choice_ind, crop_ids):  # if the current land cover is a crop
        profit_last = profits_current  # last years accumulated profit in this location
        profit_expected = (profit_signals.reshape(num_crops, 1)/yield_avg) *  location_yield # next years anticipated profit
        
    else:
        profit_last = 0
        profit_expected = np.zeros((num_crops, 1))

    return profit_last, profit_expected

def assess_profit_neighbor(j, k, crop_id_current, profits_actual_current, close_bar, crop_id_potential):
    
    radius = int(close_bar * 10)
    search_id_bound = neighbors(crop_id_current, j, k, radius)
    search_profit_bound = neighbors(profits_actual_current, j, k, radius)

    arr = np.where(search_id_bound == crop_id_potential)
    loc_arr = np.array((arr[0], arr[1])).transpose()
    
    neighbor_profit = []
    for i in range(loc_arr.shape[0]):
        neighbor_profit.append(search_profit_bound[loc_arr[i, 0], loc_arr[i, 1]])

    if len(neighbor_profit) == 0:
        neighbor_profit_avg = 0
    else:
        neighbor_profit_avg = sum(neighbor_profit) / len(neighbor_profit)

    
    return neighbor_profit_avg

def define_familiarity(d2egy_crop, switch, close_bar, far_bar):

    if 0 <= d2egy_crop < close_bar:  # close to miscanthus farmer, switching easier
        k_f = 1             # k for familiarity of miscanthus

    elif (d2egy_crop <= far_bar) & (d2egy_crop >= close_bar):  # middle to miscanthus farmer, switching neutral
        k_f = 2

    else:  # far to miscanthus farmer, switching notlikely
        k_f = 0

    alpha_f = switch[k_f][0]
    beta_f = switch[k_f][1]

    return alpha_f, beta_f

def profit_maximizer(SizeAlpha, SizeBeta, FamiliarAlpha2bg, FamiliarBeta2bg, FamiliarAlpha2chp, FamiliarBeta2chp, fmin, fmax, n, profits_current, vec_crops, vec_profit_p, profit_neighborbg, profit_neighborchp):
    # Decide which crop and associated profit to pick out of N options.

    AccRej = np.zeros(vec_crops.shape, dtype='int')

    # only consider the adoption decision for bioenergy crops, farmer will not switch back to fodder crops
    AccRej[2] = decide2switch(SizeAlpha, SizeBeta, FamiliarAlpha2bg, FamiliarBeta2bg, 
            fmin, fmax, n, profits_current, vec_profit_p[2], profit_neighborbg, seed=False)


    AccRej[3] = decide2switch(SizeAlpha, SizeBeta, FamiliarAlpha2chp, FamiliarBeta2chp, 
            fmin, fmax, n, profits_current, vec_profit_p[3], profit_neighborchp, seed=False)


    # Find the Crop IDs and associated profits that were returned as "viable": decide2switch came back as "yes" == 1
    ViableCrops = vec_crops[AccRej == 1]
    ViableProfits = vec_profit_p[AccRej == 1]

    if (ViableCrops.size == 0):
        return -1, -1

    # Find the maximum anticipated profit and the crop IDs associated with that
    MaxProfit = ViableProfits.max()
    MaxProfitCrop = ViableCrops[ViableProfits == MaxProfit]


    if (MaxProfitCrop.size > 1):
        ViableCrops = MaxProfitCrop
        ViableProfits = ViableProfits[ViableProfits == MaxProfit]
        # rule = False  # Switch rule to make the algorithm using the random option

        indChoice = np.random.choice(np.arange(ViableCrops.size), size=1)
        CropChoice = ViableCrops[indChoice]
        ProfitChoice = ViableProfits[indChoice]

    else:  # Return crop with largest profit
        CropChoice = MaxProfitCrop
        ProfitChoice = MaxProfit

    # Return the crop choice and associated profit
    return CropChoice, ProfitChoice


def make_choice(crop_id_last, profit_last, crop_choice, profit_choice, success_rate, fail_cost, seed=False):
    """ Compare the crop choice with associated profit, set the new crop ID if switching, add variability to the
            anticipated profit
    """

    if seed:

        try:
            seed_val
        except NameError:
            print("Random seed needs to be initialized using the CropDecider.DefineSeed() Function")

        np.random.seed(seed_val)

    # Check if return values indicate the farmer shouldn't switch
    if (crop_choice == -1) and (profit_choice == -1):
        crop_id_next = crop_id_last
        # profit_act = profit_last + np.random.normal(loc=0.0, scale=1000.0, size=(1, 1, 1))  # this years actual profit
        profit_act = profit_last

        # print(crop_id_last)

    else:  # switch to the new crop and add variability to resulting profit

        # successful or fail, if fail, profit is negative = initial cost
        rand_val = np.random.rand(1)

        if (rand_val < success_rate):  # successful

            crop_id_next = crop_choice
            # profit_act = profit_choice + np.random.normal(loc=0.0, scale=1000.0, size=(1, 1, 1))
            profit_act = profit_choice

        else: # failed

            crop_id_next = crop_choice

            profit_act = fail_cost

    return crop_id_next, profit_act
