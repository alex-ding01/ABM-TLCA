import numpy as np
import pandas as pd
import geopandas as gpd
import yaml

import crop_functions.crop_decider as crpdec


class ConfigReader:

    # keys found in the configuration file
    NR = 'nr'
    Y_AVERAGE = 'y_average'

    F_YD_FILE = 'f_yd_file'


    F_INIT_LC_FILE = 'f_init_lc_file'
    F_KEY_FILE = 'f_key_file'
    F_PROFITS_FILE = 'f_profits_file'
    F_DCF_FILE = 'f_dcf_file'
    F_FADN_FILE = 'f_fadn_file'
    F_LCA = 'f_lca'
    F_METADATA = 'f_metadata'
    NT = 'nt'
    NI = 'ni'
    BOUNDARY = 'boundary'
    PRODUCT = 'product'
    Y_A = 'y_a'
    Y_B = 'y_b'
    Y_C = 'y_c'
    Y_D = 'y_d'

    ENERGY_USE_IDS= 'energy_use_ids'
    C0 = 'c0'
    CA = 'ca'
    CF = 'cf'
    PT = 'pt'
    THETA = 'theta'
    ALPHA = 'alpha'    
    SUBSIDY_FARM = 'subsidy_farm'
    SUBSIDY_GREEN = 'subsidy_green'
    FAIL_COST = 'fail_cost'

    SWITCH_PARAMS_SIZE = 'switch_params_size'
    SUCCESS_RATE = 'success_rate'
    ATTR = 'attr'
    FMIN = 'fmin'
    FMAX = 'fmax'
    N = 'n'
    SCALE = 'scale'
    CROP_SEED_SIZE = 'crop_seed_size'
    TARGET_YR = 'initialization_yr'
    CLOSE_BAR = 'close_bar'
    FAR_BAR = 'far_bar'

    FARM_SIZE = 'farm_size'

    OUTPUT_DIR = 'output_directory'

    # in scenario 4-9, demonstration area are in one of the regions (default as north and south). 
    N_REGIONS_ROW = 'n_regions_row' 
    N_REGIONS_COL = 'n_regions_col'

    BATCH_LENGTH = 'batch_length'
    BATCH_WIDTH = 'batch_width'
    REGIONS = 'regions'
    N_RANDOM_BLOCKS = 'n_random_blocks'
    N_RANDOM_CELLS = 'n_random_cells'


    def __init__(self, config_file):

        c = self.read_yaml(config_file)

        self.nr = c[ConfigReader.NR]

        self.y_average = c[ConfigReader.Y_AVERAGE]

        self.f_yd_file = c[ConfigReader.F_YD_FILE]

        self.lca_folder = c[ConfigReader.F_LCA]

        self.f_metadata = c[ConfigReader.F_METADATA]

        self.f_init_lc_file = c[ConfigReader.F_INIT_LC_FILE]

        self.profits_file = c[ConfigReader.F_PROFITS_FILE]

        self.dcf_file = c[ConfigReader.F_DCF_FILE]

        self.key_file = pd.read_csv(c[ConfigReader.F_KEY_FILE])

        self.output_dir = c[ConfigReader.OUTPUT_DIR]

        # self.f_nut = c[ConfigReader.F_NUT]

        self.fadn_file = c[ConfigReader.F_FADN_FILE]

        self.Nt = c[ConfigReader.NT]

        self.Ni = c[ConfigReader.NI]

        self.boundary = c[ConfigReader.BOUNDARY]

        self.product = c[ConfigReader.PRODUCT]

        self.y_a = c[ConfigReader.Y_A]

        self.y_b = c[ConfigReader.Y_B]

        self.y_c = c[ConfigReader.Y_C]

        self.y_d = c[ConfigReader.Y_D]

        self.energy_use_ids = c[ConfigReader.ENERGY_USE_IDS]
        self.c0 = c[ConfigReader.C0]
        self.ca = c[ConfigReader.CA]
        self.cf = c[ConfigReader.CF]
        self.pt = c[ConfigReader.PT]
        self.theta = c[ConfigReader.THETA]
        self.alpha = c[ConfigReader.ALPHA]
        self.subsidy_green = c[ConfigReader.SUBSIDY_GREEN]
        self.fail_cost = c[ConfigReader.FAIL_COST]


        # set agent switching parameters (alpha, beta) [[switching averse], [switching tolerant]]
        self.switch_size = np.array(c[ConfigReader.SWITCH_PARAMS_SIZE])


        self.success_rate = c[ConfigReader.SUCCESS_RATE]

        # boolean that sets whether to base switching parameters on age and tenure (True) or not
        self.attr = c[ConfigReader.ATTR]

        # fraction of current profit at which the CDF is zero and one, and number of points to generate
        self.fmin = c[ConfigReader.FMIN]
        self.fmax = c[ConfigReader.FMAX]
        self.n = c[ConfigReader.N]

        self.scale = c[ConfigReader.SCALE]
        
        crpdec.define_seed(c[ConfigReader.CROP_SEED_SIZE])

        # target year
        self.target_year = c[ConfigReader.TARGET_YR]

        # close bar for very familiar
        self.close_bar = c[ConfigReader.CLOSE_BAR]

        # far bar for not familiar
        self.far_bar = c[ConfigReader.FAR_BAR]

        self.farm_size = c[ConfigReader.FARM_SIZE]

        self.n_regions_row = c[ConfigReader.N_REGIONS_ROW]

        self.n_regions_col = c[ConfigReader.N_REGIONS_COL]
        self.batch_length = c[ConfigReader.BATCH_LENGTH]
        self.batch_width = c[ConfigReader.BATCH_WIDTH]
        self.regions = c[ConfigReader.REGIONS]
        self.n_random_blocks = c[ConfigReader.N_RANDOM_BLOCKS]
        self.n_random_cells = c[ConfigReader.N_RANDOM_CELLS]


    @staticmethod
    def read_yaml(config_file):
        """Read the YAML config file to a dictionary.

        :param config_file:             Full path with file name and extension to the input config file.

        :return:                        YAML dictionary-like object

        """

        with open(config_file) as f:
            return yaml.safe_load(f)
