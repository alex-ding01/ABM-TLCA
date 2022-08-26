
import pandas as pd
import geopandas as gpd
from functools import reduce
import brightway2 as bw  # conda install -y -q -c conda-forge -c cmutel -c haasad brightway2
from bw2io import BW2Package
import os
import numpy as np
import rasterio
import affine
from pandarus import *
import zipfile

crs_licor = 31370

def get_initials(name):
    names = name.split(' ')
    initials = ''.join([f"{n[0]}" for n in names])
    return initials.lower()

    
def read_yieldnc(yield_fp, cropname):

    '''

    :param cropname: potential energy crop's name, read from a .nc file
    :param rcfname: regional characterization factors, read from a shapefile
    :return: yield and rCF intersected grid, present as a shapefile
    '''

    with rasterio.open(
            r'NetCDF:' + yield_fp + ':' + cropname) as src:
        yield_data = src.read(1)
    profile = {'driver': 'GTiff',
               'dtype': 'float32',
               'nodata': 1.0000000200408773e+20,
               'width': 720,
               'height': 360,
               'count': 1,
               'crs': {'init': 'epsg:4326'},
               'transform': affine.Affine(0.5, 0.0, -180.0, 0.0, -0.5, 90.0)}

    # Write raster layer as geotiff
    with rasterio.open('%s.tif' % cropname, "w", **profile) as raster:
        raster.write(yield_data, 1)

    # # show yield map as rasterfile
    # with rasterio.open('%s.tif' % cropname) as raster:
    #     img = rp.show(raster)

    # convert raster file to vector using pandarus
    # resulted data are stored in: /Users/tianran/Library/Application Support/pandarus/raster-conversion/
    crop_json = convert_to_vector('%s.tif' % cropname, band=1)
    crop_data = gpd.read_file(crop_json)
    # crs = {'init': f'epsg:{crs}'}
    crop_gpd = gpd.GeoDataFrame(crop_data, geometry=crop_data.geometry)

    return crop_gpd


def overlay3(df1, df2, df3):
    o1 = gpd.overlay(df1, df2)
    o2 = gpd.overlay(o1, df3)
    
    return o2

def ecoinvent_set_up():
    print('Brightway2 set up')
    bw.projects.set_current('tef_calculation3')
    bw.bw2setup()
    bidb_infunc = bw.Database('biosphere3')

    if 'ecoinvent 3.5 cutoff' in bw.databases:
        print('Database has already been imported')

    else:
        ei35 = bw.SingleOutputEcospold2Importer(
            "data/ecoinvent 3.5_cutoff_ecoSpold02/datasets",
            "ecoinvent 3.5 cutoff")
        ei35.apply_strategies()
        ei35.statistics()
        ei35.write_database()

    eidb_infunc = bw.Database('ecoinvent 3.5 cutoff')
    # bw.create_core_migrations()

    # Read ecoinvent database option 1
    # eidl.get_ecoinvent()
    # # Username: ULB
    # # password: polytechLCA8

    # Read ecoinvent database option 2

    return [bidb_infunc, eidb_infunc]



def emis_cf(file_path):
    """

    :param file_path: LCIA method file path, downloaded from http://www.impactworldplus.org/en/writeToFile.php
    :return: cf for each emissions for mid-point indicators
    """
    iw_mid = BW2Package.import_file(file_path)
    lca_method = iw_mid

    lcia_name = []
    emis_name = []
    emis_cate = []
    dcf = []

    for m in lca_method:
        
        for cf in m.load():  # emission index and cf
          
            emis = bw.get_activity(cf[0])  # emission dictionary information for below

            m_name = str(m).split(": ")[-1]
            e_name = emis.as_dict()['name']  # emission name
            e_cate = emis.as_dict()['categories']  # emission category
            cf = cf[1]  # emission CF

            lcia_name.append(m_name)
            emis_name.append(e_name)
            emis_cate.append(e_cate)
            dcf.append(cf)

    df = pd.DataFrame()
    df['lcia_name'] = lcia_name
    df['emis_name'] = emis_name
    df['emis_cate'] = emis_cate
    df['dcf'] = dcf

    return df


def cal_tef_misc(inv_misc_cult, iw_method):

    inpt_name = []
    emis_name = []
    emis_sorc = []
    emis_cate = []
    emis_amnt = []

    mx_territory_input_new = [
    'diesel, burned in agricultural machinery'
    ]
    of_territory_input_new = [
    'calcium ammonium nitrate production',
    'triple superphosphate production',
    'potassium chloride production',
    'pesticide production, unspecified'
    ]


    for tech in inv_misc_cult.technosphere():

        inpt_name_value = tech.input.as_dict()['name']

        inpt_amnt_value = tech['amount']  # calculate the amount for 1 hectare of product

        for bio in tech.input.biosphere():

            inpt_name_value = 'na'
            emis_name_value = bio.as_dict()['name']
            emis_cate_value = bio.input.as_dict()['categories']
            emis_amnt_value = bio.as_dict()['amount']  # calculate the amount for 1 hectare of product
            emis_sorc_value = 'in-territory emission'

            inpt_name.append(inpt_name_value)
            emis_name.append(emis_name_value)
            emis_sorc.append(emis_sorc_value)
            emis_cate.append(emis_cate_value)
            emis_amnt.append(emis_amnt_value)

            # the inputs type in-off-mix information

        if tech['name'] in of_territory_input_new:

            for m in iw_method[:1]:

                lcia_name_value = m[2]

                lca = bw.LCA({tech.input: inpt_amnt_value}, m)
                lca.lci()
                lca.lcia()
                lcia_dfsc_value = lca.score

                for e in lca.biosphere_dict:
                    emis = bw.get_activity(e)
                    emis_name_value = emis.as_dict()['name']
                    emis_cate_value = emis.as_dict()['categories']
                    emis_sorc_value = 'off-territory emission'

                    emis_row = lca.biosphere_dict[emis]
                    emis_amnt_value = lca.inventory[emis_row, :].sum()

                    inpt_name.append(inpt_name_value)
                    emis_name.append(emis_name_value)
                    emis_sorc.append(emis_sorc_value)
                    emis_cate.append(emis_cate_value)
                    emis_amnt.append(emis_amnt_value)

        else:

            for tech2 in tech.input.technosphere():

                tech_name_value = tech2.input.as_dict()['name']
                tech2_amnt = tech2['amount']


                            # calculate the LCA score for these technosphere input of the mix-input
                for m in iw_method[:1]:
                    lcia_name_value = m[2]

                    lca = bw.LCA({tech2.input: inpt_amnt_value * tech2_amnt}, m)
                    lca.lci()
                    lca.lcia()
                    lcia_dfsc_value = lca.score


                    for e in lca.biosphere_dict:
                        emis = bw.get_activity(e)
                        emis_name_value = emis.as_dict()['name']
                        emis_cate_value = emis.as_dict()['categories']
                        emis_sorc_value = 'off-territory emission'

                        emis_row = lca.biosphere_dict[emis]
                        emis_amnt_value = lca.inventory[emis_row, :].sum()

                        inpt_name.append(inpt_name_value)
                        emis_name.append(emis_name_value)
                        emis_sorc.append(emis_sorc_value)
                        emis_cate.append(emis_cate_value)
                        emis_amnt.append(emis_amnt_value)

                for bio2 in tech2.input.biosphere():
                    emis_name_value = bio2.as_dict()['name']
                    emis_sorc_value = 'in-territory emission'
                    emis_cate_value = bio2.input.as_dict()['categories']
                    emis_amnt_value = inpt_amnt_value * tech2_amnt * bio2['amount']

                    inpt_name.append(inpt_name_value)
                    emis_name.append(emis_name_value)
                    emis_sorc.append(emis_sorc_value)
                    emis_cate.append(emis_cate_value)
                    emis_amnt.append(emis_amnt_value)

    tef_new = pd.DataFrame()
    tef_new['inpt_name'] = inpt_name
    tef_new['emis_name'] = emis_name
    tef_new['emis_sorc'] = emis_sorc
    tef_new['emis_cate'] = emis_cate
    tef_new['emis_amnt'] = emis_amnt

    return tef_new


def cal_tef_pre(inv, iw_method):

    inpt_name = []
    emis_name = []
    emis_sorc = []
    emis_cate = []
    emis_amnt = []

    in_territory_input_new = [
    'maize seed production, for sowing',
    'poultry manure, fresh, Recycled Content cut-off',
    'manure, liquid, swine, Recycled Content cut-off',
    'grass seed production, Swiss integrated production, for sowing'
    ]
    mx_territory_input_new = [
    'transport, tractor and trailer, agricultural',
    'tillage, harrowing, by rotary harrow',
    'tillage, ploughing',
    'application of plant protection product, by field sprayer',
    'chopping, maize',
    'fertilising, by broadcaster',
    'sowing',
    'fodder loading, by self-loading trailer',
    'hoeing',
    'operation, housing system, cattle, loose',
    'operation, housing system, cattle, tied',
    'operation, housing system, pig, fully-slatted floor',
    'Bale loading',
    'Baling',
    'Haying, by rotary tedder',
    'Mowing, by rotary mower', 
    'Swath, by rotary windrower', 
    'tillage, cultivating, chiselling'
    ]

    for tech in inv.technosphere():
        inpt_name_value = tech.input.as_dict()['name']
        inpt_amnt_value = tech['amount']  # calculate the amount for 1 hectare of product

        for bio2 in tech.input.biosphere():
            emis_name_value = bio2.as_dict()['name']
            emis_sorc_value = 'in-territory emission'
            emis_cate_value = bio2.input.as_dict()['categories']
            emis_amnt_value = inpt_amnt_value * bio2.as_dict()['amount']

            inpt_name.append('na')
            emis_name.append(emis_name_value)
            emis_sorc.append(emis_sorc_value)
            emis_cate.append(emis_cate_value)
            emis_amnt.append(emis_amnt_value)
    
        for tech2 in tech.input.technosphere():
            tech2_name_value = tech2['name']
            tech2_amnt = tech2['amount']
            
            if tech2['name'] in in_territory_input_new:

                for m in iw_method[:1]:
                    
                    lcia_name_value = m[2]

                    lca = bw.LCA({tech2.input: inpt_amnt_value * tech2_amnt}, m)
                    lca.lci()
                    lca.lcia()
                    lcia_dfsc_value = lca.score

                    for e in lca.biosphere_dict:
                        emis = bw.get_activity(e)
                        emis_name_value = emis.as_dict()['name']
                        emis_cate_value = emis.as_dict()['categories']
                        emis_sorc_value = 'in-territory emission'

                        emis_row = lca.biosphere_dict[emis]
                        emis_amnt_value = lca.inventory[emis_row, :].sum()
                        
                        inpt_name.append(inpt_name_value)
                        emis_name.append(emis_name_value)
                        emis_sorc.append(emis_sorc_value)
                        emis_cate.append(emis_cate_value)
                        emis_amnt.append(emis_amnt_value)
                        
            elif tech2['name'] in mx_territory_input_new:
                
                for tech3 in tech2.input.technosphere():
                    tech3_amnt = tech3['amount']
                    
                    
                    for m in iw_method[:1]:
                    
                        lcia_name_value = m[2]

                        lca = bw.LCA({tech3.input: inpt_amnt_value * tech2_amnt * tech3_amnt}, m)
                        lca.lci()
                        lca.lcia()
                        lcia_dfsc_value = lca.score

                        for e in lca.biosphere_dict:
                            emis = bw.get_activity(e)
                            emis_name_value = emis.as_dict()['name']
                            emis_cate_value = emis.as_dict()['categories']
                            emis_sorc_value = 'off-territory emission'

                            emis_row = lca.biosphere_dict[emis]
                            emis_amnt_value = lca.inventory[emis_row, :].sum()
                            
                            inpt_name.append(inpt_name_value)
                            emis_name.append(emis_name_value)
                            emis_sorc.append(emis_sorc_value)
                            emis_cate.append(emis_cate_value)
                            emis_amnt.append(emis_amnt_value)
                for bio3 in tech2.input.biosphere():
                        emis_name_value = bio3.as_dict()['name']
                        emis_sorc_value = 'in-territory emission'
                        emis_cate_value = bio3.input.as_dict()['categories']
                        emis_amnt_value = inpt_amnt_value * tech2_amnt * bio3.as_dict()['amount']

                        
                        inpt_name.append(inpt_name_value)
                        emis_name.append(emis_name_value)
                        emis_sorc.append(emis_sorc_value)
                        emis_cate.append(emis_cate_value)
                        emis_amnt.append(emis_amnt_value)
                        
            else:
                for m in iw_method[:1]:
                    
                    lcia_name_value = m[2]

                    lca = bw.LCA({tech2.input: inpt_amnt_value * tech2_amnt}, m)
                    lca.lci()
                    lca.lcia()
                    lcia_dfsc_value = lca.score

                    for e in lca.biosphere_dict:
                        emis = bw.get_activity(e)
                        emis_name_value = emis.as_dict()['name']
                        emis_cate_value = emis.as_dict()['categories']
                        emis_sorc_value = 'off-territory emission'

                        emis_row = lca.biosphere_dict[emis]
                        emis_amnt_value = lca.inventory[emis_row, :].sum()

                        inpt_name.append(inpt_name_value)
                        emis_name.append(emis_name_value)
                        emis_sorc.append(emis_sorc_value)
                        emis_cate.append(emis_cate_value)
                        emis_amnt.append(emis_amnt_value)
            
    tef_new = pd.DataFrame()
    tef_new['inpt_name'] = inpt_name
    tef_new['emis_name'] = emis_name
    tef_new['emis_sorc'] = emis_sorc
    tef_new['emis_cate'] = emis_cate
    tef_new['emis_amnt'] = emis_amnt

    return tef_new

# def cal_tef_gras(inv_gras, iw_method):
#     pass






def cal_tef1(crops, iw_method):
    crop_name = []
    inpt_name = []
    inpt_type = []
    inpt_amnt = []
    emis_name = []
    emis_sorc = []
    emis_cate = []
    emis_amnt = []

    in_territory_input_new = [
        'beet seed, Swiss integrated production, for sowing',
        'miscanthus rhizome, for planting',
    ]
    mx_territory_input_new = [
        'application of plant protection product, by field sprayer',
        'combine harvesting',
        'fertilising, by broadcaster',
        'fodder loading, by self-loading trailer',
        'hoeing',
        'mulching',
        'potato planting',
        'sowing',
        'tillage, cultivating, chiselling',
        'tillage, currying, by weeder',
        'tillage, harrowing, by rotary harrow',
        'tillage, harrowing, by spring tine harrow',
        'tillage, ploughing',
        'tillage, rolling',
        'tillage, rotary cultivator',
        'transport, tractor and trailer, agricultural',
        'chopping, maize',
        'tillage, ploughing',
        'drying of bread grain, seed and legumes',
        'irrigation'
    ]
    of_territory_input_new = [
        'glyphosate',
        'nitrogen fertiliser, as N',
        'phosphate rock, as P2O5, beneficiated, dry',
        'potassium sulfate, as K2O',
        'nitrogen fertiliser, as N',
        'acetamide-anillide-compound, unspecified',
        'pesticide, unspecified',
        'lime',
        'building, hall, wood construction',
        'ammonium nitrate, as N',
        'cyclic N-compound',
        'diphenylether-compound',
        'nitrogen fertiliser, as N',
        'pesticide, unspecified',
        'phosphate fertiliser, as P2O5',
        'potassium chloride, as K2O',
        'pyrethroid-compound',
        'magnesium oxide'
    ]

    test_in_te_new = []
    test_of_te_new = []
    test_mx_te_new = []

    for crop in crops:
        crop_name_value = crop['name']
        print()
        print(crop_name_value)

        # three energy crops have different input name
        # direct emissions
        for bio in crop.biosphere():
            inpt_name_value = 'na'
            inpt_type_value = 'direct emission'
            inpt_amnt_value = 'na'

            emis_name_value = bio.as_dict()['name']
            emis_cate_value = bio.input.as_dict()['categories']
            emis_amnt_value = bio.as_dict()['amount']  # calculate the amount for 1 hectare of product
            emis_sorc_value = 'in-territory emission'

            crop_name.append(crop_name_value)
            inpt_name.append(inpt_name_value)
            inpt_type.append(inpt_type_value)
            inpt_amnt.append(inpt_amnt_value)

            emis_name.append(emis_name_value)
            emis_sorc.append(emis_sorc_value)
            emis_cate.append(emis_cate_value)
            emis_amnt.append(emis_amnt_value)

            # technosphere
        for tech in crop.technosphere():

            inpt_name_value = tech.input.as_dict()['name']

            inpt_amnt_value = tech['amount']  # calculate the amount for 1 hectare of product

            # the inputs type in-off-mix information
            if tech['name'] in in_territory_input_new:
                test_in_te_new.append(tech['name'])
                inpt_type_value = 'in-territory input'

                for m in iw_method[:1]:

                    lcia_name_value = m[2]

                    lca = bw.LCA({tech.input: inpt_amnt_value}, m)
                    lca.lci()
                    lca.lcia()
                    lcia_dfsc_value = lca.score

                    for e in lca.biosphere_dict:
                        emis = bw.get_activity(e)
                        emis_name_value = emis.as_dict()['name']
                        emis_cate_value = emis.as_dict()['categories']
                        emis_sorc_value = 'in-territory emission'

                        emis_row = lca.biosphere_dict[emis]
                        emis_amnt_value = lca.inventory[emis_row, :].sum()

                        crop_name.append(crop_name_value)
                        inpt_name.append(inpt_name_value)
                        inpt_type.append(inpt_type_value)
                        inpt_amnt.append(inpt_amnt_value)

                        emis_name.append(emis_name_value)
                        emis_sorc.append(emis_sorc_value)
                        emis_cate.append(emis_cate_value)
                        emis_amnt.append(emis_amnt_value)

            elif tech['name'] in of_territory_input_new:
                test_of_te_new.append(tech['name'])
                inpt_type_value = 'off-territory input'
                for m in iw_method[:1]:

                    lcia_name_value = m[2]

                    lca = bw.LCA({tech.input: inpt_amnt_value}, m)
                    lca.lci()
                    lca.lcia()
                    lcia_dfsc_value = lca.score

                    for e in lca.biosphere_dict:
                        emis = bw.get_activity(e)
                        emis_name_value = emis.as_dict()['name']
                        emis_cate_value = emis.as_dict()['categories']
                        emis_sorc_value = 'off-territory emission'

                        emis_row = lca.biosphere_dict[emis]
                        emis_amnt_value = lca.inventory[emis_row, :].sum()

                        crop_name.append(crop_name_value)
                        inpt_name.append(inpt_name_value)
                        inpt_type.append(inpt_type_value)
                        inpt_amnt.append(inpt_amnt_value)
                        emis_name.append(emis_name_value)
                        emis_sorc.append(emis_sorc_value)
                        emis_cate.append(emis_cate_value)
                        emis_amnt.append(emis_amnt_value)

            else:
                test_mx_te_new.append(tech['name'])
                inpt_type_value = 'mix-territory input'

                for tech2 in tech.input.technosphere():

                    tech_name_value = tech2.input.as_dict()['name']
                    tech2_amnt = tech2['amount']

                    for tech3 in tech2.input.technosphere():

                        tech3_amnt = tech3['amount']

                        tect_amnt = tech2_amnt * tech3_amnt

                        # calculate the LCA score for these technosphere input of the mix-input
                        for m in iw_method[:1]:
                            lcia_name_value = m[2]

                            lca = bw.LCA({tech3.input: inpt_amnt_value * tect_amnt}, m)
                            lca.lci()
                            lca.lcia()
                            lcia_dfsc_value = lca.score

                            for e in lca.biosphere_dict:
                                emis = bw.get_activity(e)
                                emis_name_value = emis.as_dict()['name']
                                emis_cate_value = emis.as_dict()['categories']
                                emis_sorc_value = 'off-territory emission'

                                emis_row = lca.biosphere_dict[emis]
                                emis_amnt_value = lca.inventory[emis_row, :].sum()

                                crop_name.append(crop_name_value)
                                inpt_name.append(inpt_name_value)
                                inpt_type.append(inpt_type_value)
                                inpt_amnt.append(inpt_amnt_value)
                                emis_name.append(emis_name_value)
                                emis_sorc.append(emis_sorc_value)
                                emis_cate.append(emis_cate_value)
                                emis_amnt.append(emis_amnt_value)

                    for bio3 in tech2.input.biosphere():
                        emis_name_value = bio3.as_dict()['name']
                        emis_sorc_value = 'in-territory emission'
                        emis_cate_value = bio3.input.as_dict()['categories']
                        emis_amnt_value = inpt_amnt_value * tech2_amnt * bio3['amount']

                        crop_name.append(crop_name_value)
                        inpt_name.append(inpt_name_value)
                        inpt_type.append(inpt_type_value)
                        inpt_amnt.append(inpt_amnt_value)
                        emis_name.append(emis_name_value)
                        emis_sorc.append(emis_sorc_value)
                        emis_cate.append(emis_cate_value)
                        emis_amnt.append(emis_amnt_value)

    tef_new = pd.DataFrame()
    tef_new['crop_name'] = crop_name
    tef_new['inpt_name'] = inpt_name
    tef_new['inpt_type'] = inpt_type
    tef_new['inpt_amnt'] = inpt_amnt
    tef_new['emis_name'] = emis_name
    tef_new['emis_sorc'] = emis_sorc
    tef_new['emis_cate'] = emis_cate
    tef_new['emis_amnt'] = emis_amnt  # emission amount for each input/ha

    # tef_new.to_csv("/Users/tianran/Desktop/janustd/data_td/Input4stlca/misc_emis_pkg.csv")

    return tef_new


def cal_tef2(crop):

    # only for maize

    in_territory_input = [
        'green manure, Swiss integrated production, until January',
        'green manure, Swiss integrated production, until February',
        'green manure, Swiss integrated production, until April',
        'maize seed, Swiss integrated production, for sowing',
        'pea seed, for sowing',
        'fava bean seed, for sowing',
        'fodder beet seed, for sowing']
    of_territory_input = [
        '[thio]carbamate-compound',
        'ammonium nitrate, as N',
        'phosphate fertiliser, as P2O5',
        'pesticide, unspecified',
        'metolachlor',
        'potassium sulfate, as K2O',
        'atrazine',
        'potassium chloride, as K2O',
        'nitrogen fertiliser, as N',
        'urea, as N',
        'ammonium sulfate, as N',
        'phosphate rock, as P2O5, beneficiated, dry',
        'glyphosate',
        'ammonium nitrate, as N',
        'benzo[thia]diazole-compound',
        'dinitroaniline-compound',
        'pesticide, unspecified',
        'phosphate rock, as P2O5, beneficiated, dry',
        'nitrile-compound',
        'diphenylether-compound',
        'triazine-compound, unspecified',
        'organophosphorus-compound, unspecified',
        'ammonium sulfate, as N',
        'benzimidazole-compound'
    ]
    mx_territory_input_market = [
        'irrigation',
        'transport, tractor and trailer, agricultural'
    ]

    crop_name = []
    crop_yild = []

    inpt_name = []
    inpt_type = []  # in-territory input, mix-territory input, off-territory input
    inpt_amnt = []

    emis_name = []
    emis_sorc = []  # define if it is in-territory emission or off-territory emission
    #                 'in-territory emission', 'off-territory emission'
    emis_cate = []
    emis_amnt = []

    test_in_te = []
    test_of_te = []
    test_mx_te = []
    test_mx_te_market = []  # the input for mix input is from market

    yild = yield_default[2]
    crop_name_value = 'maize silage production, Swiss integrated production, intensive'

    # direct emissions
    for bio in crop.biosphere():
        inpt_name_value = 'na'
        inpt_type_value = 'direct emission'
        inpt_amnt_value = 'na'

        emis_name_value = bio.as_dict()['name']
        emis_cate_value = bio.input.as_dict()['categories']
        emis_amnt_value = bio.as_dict()['amount'] * yild  # calculate the amount for 1 hectare of product
        emis_sorc_value = 'in-territory emission'

        crop_name.append(crop_name_value)
        crop_yild.append(yild)
        inpt_name.append(inpt_name_value)
        inpt_type.append(inpt_type_value)
        inpt_amnt.append(inpt_amnt_value)

        emis_name.append(emis_name_value)
        emis_sorc.append(emis_sorc_value)
        emis_cate.append(emis_cate_value)
        emis_amnt.append(emis_amnt_value)

    # technosphere
    for tech in crop.technosphere():

        inpt_name_value = tech.input.as_dict()['name']

        inpt_amnt_value = tech['amount'] * yild  # calculate the amount for 1 hectare of product

        # the inputs type in-off-mix information
        if tech['name'] in in_territory_input:
            test_in_te.append(tech['name'])
            inpt_type_value = 'in-territory input'
            for m in iw_method[:1]:

                lcia_name_value = m[2]
                # the input name & amount

                lca = bw.LCA({tech.input: inpt_amnt_value}, m)
                lca.lci()
                lca.lcia()
                lcia_dfsc_value = lca.score

                for e in lca.biosphere_dict:
                    emis = bw.get_activity(e)
                    emis_name_value = emis.as_dict()['name']
                    emis_cate_value = emis.as_dict()['categories']
                    emis_sorc_value = 'in-territory emission'

                    emis_row = lca.biosphere_dict[emis]
                    emis_amnt_value = lca.inventory[emis_row, :].sum()

                    crop_name.append(crop_name_value)
                    crop_yild.append(yild)
                    inpt_name.append(inpt_name_value)
                    inpt_type.append(inpt_type_value)
                    inpt_amnt.append(inpt_amnt_value)

                    emis_name.append(emis_name_value)
                    emis_sorc.append(emis_sorc_value)
                    emis_cate.append(emis_cate_value)
                    emis_amnt.append(emis_amnt_value)

        elif tech['name'] in of_territory_input:
            test_of_te.append(tech['name'])
            inpt_type_value = 'off-territory input'
            for m in iw_method[:1]:

                lcia_name_value = m[2]

                lca = bw.LCA({tech.input: inpt_amnt_value}, m)
                lca.lci()
                lca.lcia()
                lcia_dfsc_value = lca.score

                for e in lca.biosphere_dict:
                    emis = bw.get_activity(e)
                    emis_name_value = emis.as_dict()['name']
                    emis_cate_value = emis.as_dict()['categories']
                    emis_sorc_value = 'off-territory emission'

                    emis_row = lca.biosphere_dict[emis]
                    emis_amnt_value = lca.inventory[emis_row, :].sum()

                    crop_name.append(crop_name_value)
                    crop_yild.append(yild)
                    inpt_name.append(inpt_name_value)
                    inpt_type.append(inpt_type_value)
                    inpt_amnt.append(inpt_amnt_value)
                    emis_name.append(emis_name_value)
                    emis_sorc.append(emis_sorc_value)
                    emis_cate.append(emis_cate_value)
                    emis_amnt.append(emis_amnt_value)

        elif tech['name'] in mx_territory_input_market:

            inpt_type_value = 'mix-territory input'

            test_mx_te_market.append(tech['name'])

            for tech2 in tech.input.technosphere():

                tech_name_value = tech2.input.as_dict()['name']
                tech2_amnt = tech2['amount']

                for tech3 in tech2.input.technosphere():

                    tech3_amnt = tech3['amount']

                    tect_amnt = tech2_amnt * tech3_amnt

                    # calculate the LCA score for these technosphere input of the mix-input
                    for m in iw_method[:1]:
                        lcia_name_value = m[2]

                        lca = bw.LCA({tech3.input: inpt_amnt_value * tect_amnt}, m)
                        lca.lci()
                        lca.lcia()
                        lcia_dfsc_value = lca.score

                        for e in lca.biosphere_dict:
                            emis = bw.get_activity(e)
                            emis_name_value = emis.as_dict()['name']
                            emis_cate_value = emis.as_dict()['categories']
                            emis_sorc_value = 'off-territory emission'

                            emis_row = lca.biosphere_dict[emis]
                            emis_amnt_value = lca.inventory[emis_row, :].sum()

                            crop_name.append(crop_name_value)
                            crop_yild.append(yild)
                            inpt_name.append(inpt_name_value)
                            inpt_type.append(inpt_type_value)
                            inpt_amnt.append(inpt_amnt_value)
                            emis_name.append(emis_name_value)
                            emis_sorc.append(emis_sorc_value)
                            emis_cate.append(emis_cate_value)
                            emis_amnt.append(emis_amnt_value)

                for bio3 in tech2.input.biosphere():
                    emis_name_value = bio3.as_dict()['name']
                    emis_sorc_value = 'in-territory emission'
                    emis_cate_value = bio3.input.as_dict()['categories']
                    emis_amnt_value = inpt_amnt_value * tech2_amnt * bio3['amount']

                    crop_name.append(crop_name_value)
                    crop_yild.append(yild)
                    inpt_name.append(inpt_name_value)
                    inpt_type.append(inpt_type_value)
                    inpt_amnt.append(inpt_amnt_value)
                    emis_name.append(emis_name_value)
                    emis_sorc.append(emis_sorc_value)
                    emis_cate.append(emis_cate_value)
                    emis_amnt.append(emis_amnt_value)

        else:
            test_mx_te.append(tech['name'])

            inpt_type_value = 'mix-territory input'

            for tech2 in tech.input.technosphere():  # impacts belongs to off-territory emissions

                tech2_amnt = tech2['amount']  # input amount for 1 unit process of transport

                tech_amnt_value = tech2_amnt
                # calculate the LCA score for these technosphere input of the mix-input
                for m in iw_method[:1]:
                    lcia_name_value = m[2]

                    lca = bw.LCA({tech2.input: inpt_amnt_value * tech2_amnt}, m)
                    lca.lci()
                    lca.lcia()
                    lcia_dfsc_value = lca.score

                    for e in lca.biosphere_dict:
                        emis = bw.get_activity(e)
                        emis_name_value = emis.as_dict()['name']
                        emis_cate_value = emis.as_dict()['categories']
                        emis_sorc_value = 'off-territory emission'

                        emis_row = lca.biosphere_dict[emis]
                        emis_amnt_value = lca.inventory[emis_row, :].sum()

                        crop_name.append(crop_name_value)
                        crop_yild.append(yild)
                        inpt_name.append(inpt_name_value)
                        inpt_type.append(inpt_type_value)
                        inpt_amnt.append(inpt_amnt_value)
                        emis_name.append(emis_name_value)
                        emis_sorc.append(emis_sorc_value)
                        emis_cate.append(emis_cate_value)
                        emis_amnt.append(emis_amnt_value)

            for bio2 in tech.input.biosphere():
                emis_name_value = bio2.as_dict()['name']
                emis_sorc_value = 'in-territory emission'
                emis_cate_value = bio2.input.as_dict()['categories']
                emis_amnt_value = bio2['amount'] * inpt_amnt_value

                crop_name.append(crop_name_value)
                crop_yild.append(yild)
                inpt_name.append(inpt_name_value)
                inpt_type.append(inpt_type_value)
                inpt_amnt.append(inpt_amnt_value)
                emis_name.append(emis_name_value)
                emis_sorc.append(emis_sorc_value)
                emis_cate.append(emis_cate_value)
                emis_amnt.append(emis_amnt_value)


    tef_new = pd.DataFrame()
    tef_new['crop_name'] = crop_name
    tef_new['crop_yild'] = crop_yild
    tef_new['inpt_name'] = inpt_name
    tef_new['inpt_type'] = inpt_type
    tef_new['inpt_amnt'] = inpt_amnt
    tef_new['emis_name'] = emis_name
    tef_new['emis_sorc'] = emis_sorc
    tef_new['emis_cate'] = emis_cate
    tef_new['emis_amnt'] = emis_amnt  # emission amount for each input/ha

    return tef_new


def read_tef():
    zf = zipfile.ZipFile('data/tef.zip')  # having First.csv zipped file.
    tef = pd.read_csv(zf.open('tef.csv'))

    return tef


def current_impact(name, tef_cur, cf_inds):
    
    tef = tef_cur[tef_cur.lu == name]
    tef = tef.groupby(['emis_name', 'emis_category'])[['emis_in_lu', 'emis_off_lu']].mean().reset_index()

    tef_lcia = tef.merge(cf_inds, on=['emis_name'])
    
    
    tef_lcia['lcia_score_in'] = tef_lcia.emis_in_lu * tef_lcia.dcf  # per ha
    tef_lcia['lcia_score_of'] = tef_lcia.emis_off_lu * tef_lcia.dcf  # per ha
    
    res = tef_lcia.groupby(['lcia_name'])[['lcia_score_in', 'lcia_score_of']].sum().reset_index()
    
    res['CULT_NOM'] = name

    return res

def rcf_no_nut(indicator):  # this function return regionalized CF in Walloon region

    # Read rCF (regionalised characterisation factors)
    # file paths of rCF downloaded from Impactworld+

    cf_AcidFW_fp = 'data/rCF/SHP/AcidFW_Damage_native.shp'
    cf_AcidTerr_fp = 'data/rCF/SHP/AcidTerr_Damage_native.shp'
    cf_EutroFW_fp = 'data/rCF/SHP/EutroFW_Damage_native.shp'
    cf_EutroMar_fp = 'data/rCF/SHP/EutroMar_Damage_native.shp'
    cf_LandOcc_fp = 'data/rCF/SHP/LandOcc_Damage_native.shp'
    cf_LandTrans_fp = 'data/rCF/SHP/LandTrans_Damage_native.shp'
    cf_WaterAvailab_HH_fp = 'data/rCF/SHP/WaterAvailab_HH_Damage_native.shp'

    if indicator == 'Freshwater acidification':
        cf_AcidFW = gpd.read_file(cf_AcidFW_fp)
        df = cf_AcidFW.loc[:, ['HNO3', 'NH3', 'NOX', 'SO2', 'SO4', 'geometry']]

    if indicator == 'Terrestrial acidification':
        cf_AcidTerr = gpd.read_file(cf_AcidTerr_fp)
        df = cf_AcidTerr.loc[:, ['NH3', 'NOX', 'SO2', 'geometry']]

    if indicator == 'Freshwater eutrophication':
        cf_EutroFW = gpd.read_file(cf_EutroFW_fp)
        df = cf_EutroFW.loc[:, ['CELL_ID', 'BOD', 'COD', 'PHOSPHATE', 'PHOSACID', 'PHOSPHORUS', 'PHOSPENTOX', 'geometry']]

    if indicator == 'Marine eutrophication':
        cf_EutroMar = gpd.read_file(cf_EutroMar_fp)
        df = cf_EutroMar.loc[:, ['HNO3', 'NH3', 'NOX', 'geometry']]

    if indicator == 'Land occupation, biodiversity':
        cf_LandOcc = gpd.read_file(cf_LandOcc_fp)
        df = cf_LandOcc.loc[:, ['IMPCAT', 'UNIT', 'RESOLUTION', 'ANNUALCROP', 'PASTURE', 'geometry']]

    elif indicator == 'Water availability, human health':
        cf_WaterAvailab_HH = gpd.read_file(cf_WaterAvailab_HH_fp)
        cf_WaterAvailab_HH_BE = cf_WaterAvailab_HH[cf_WaterAvailab_HH.ID_CTRY == 'BE']  # cut to belgium first
        df = cf_WaterAvailab_HH_BE

    
    nuts_fp = 'data/wl_n1/l1_wl.shp'
    nuts = gpd.read_file(nuts_fp).to_crs(epsg=4326).drop(columns=['NUTS_ID'])

    nuts_indicator = gpd.overlay(nuts, df, how='intersection').to_crs(epsg=31370)

    return nuts_indicator


def cal_fweutro_rlca(tef_df, rcf_df, lu):

    tef_df = tef_df[tef_df.lu == lu]
    tef_df = tef_df.groupby(['emis_name', 'emis_category'])[['emis_in_lu', 'emis_off_lu']].mean().reset_index()

    BOD = 'BOD5, Biological Oxygen Demand'
    COD = 'COD, Chemical Oxygen Demand'
    Phosphate = 'Phosphate'
    Pacid = 'Phosphoric acid'
    Phosphorus = 'Phosphorus'

    rcf_df['fweutro'] = rcf_df.BOD * (tef_df[tef_df.emis_name == BOD].emis_in_lu.sum()) + \
                       rcf_df.COD * (tef_df[tef_df.emis_name == COD].emis_in_lu.sum()) + \
                       rcf_df.PHOSPHATE * (tef_df[tef_df.emis_name == Phosphate].emis_in_lu.sum()) + \
                       rcf_df.PHOSACID * (tef_df[tef_df.emis_name == Pacid].emis_in_lu.sum()) + \
                       rcf_df.PHOSPHORUS * (tef_df[tef_df.emis_name == Phosphorus].emis_in_lu.sum()) 

    result = rcf_df.drop(columns=['BOD', 'COD', 'PHOSPHATE', 'PHOSACID', 'PHOSPHORUS', 'PHOSPENTOX'])
    result['lu'] = lu

    return result

def cal_fwacid_rlca(tef_df, rcf_df, lu):

    tef_df = tef_df[tef_df.lu == lu]
    tef_df = tef_df.groupby(['emis_name', 'emis_category'])[['emis_in_lu', 'emis_off_lu']].mean().reset_index()

    HNO3 = 'Nitrate'
    NH3 = 'Ammonia'
    NOX = 'Nitrogen oxides'
    SO2 = 'Sulfur dioxide'
    SO4 = 'Sulfate'

    rcf_df['fwacid'] = rcf_df.HNO3 * (tef_df[tef_df.emis_name == HNO3].emis_in_lu.sum()) + \
                       rcf_df.NH3 * (tef_df[tef_df.emis_name == NH3].emis_in_lu.sum()) + \
                       rcf_df.NOX * (tef_df[tef_df.emis_name == NOX].emis_in_lu.sum()) + \
                       rcf_df.SO2 * (tef_df[tef_df.emis_name == SO2].emis_in_lu.sum()) + \
                       rcf_df.SO4 * (tef_df[tef_df.emis_name == SO4].emis_in_lu.sum()) 

    result = rcf_df.drop(columns=['HNO3', 'NH3', 'NOX', 'SO2', 'SO4'])
    result['lu'] = lu

    return result

def cal_teracid_rlca(tef_df, rcf_df, lu):

    tef_df = tef_df[tef_df.lu == lu]
    tef_df = tef_df.groupby(['emis_name', 'emis_category'])[['emis_in_lu', 'emis_off_lu']].mean().reset_index()

    NH3 = 'Ammonia'
    NOX = 'Nitrogen oxides'
    SO2 = 'Sulfur dioxide'

    rcf_df['teracid'] =  rcf_df.NH3 * (tef_df[tef_df.emis_name == NH3].emis_in_lu.sum()) + \
                       rcf_df.NOX * (tef_df[tef_df.emis_name == NOX].emis_in_lu.sum()) + \
                       rcf_df.SO2 * (tef_df[tef_df.emis_name == SO2].emis_in_lu.sum()) 

    result = rcf_df.drop(columns=['NH3', 'NOX', 'SO2'])
    result['lu'] = lu

    return result

def cal_maeutro_rlca(tef_df, rcf_df, lu):

    tef_df = tef_df[tef_df.lu == lu]
    tef_df = tef_df.groupby(['emis_name', 'emis_category'])[['emis_in_lu', 'emis_off_lu']].mean().reset_index()

    NH3 = 'Ammonia'
    NOX = 'Nitrogen oxides'

    rcf_df['maeutro'] =  rcf_df.NH3 * (tef_df[tef_df.emis_name == NH3].emis_in_lu.sum()) + \
                       rcf_df.NOX * (tef_df[tef_df.emis_name == NOX].emis_in_lu.sum())

    result = rcf_df.drop(columns=['NH3', 'NOX'])
    result['lu'] = lu

    return result

def cal_tef3(process_name, iw_method):
    proc_name = []
    inpt_name = []
    inpt_type = []
    inpt_amnt = []
    emis_name = []
    emis_sorc = []
    emis_cate = []
    emis_amnt = []

    in_territory_input_new = [
        'electricity, from municipal waste incineration to generic market for electricity, medium voltage',
        'heat, from municipal waste incineration to generic market for heat district or industrial, other than natural gas'
    ]
    mx_territory_input_new = [
        'transport, tractor and trailer, agricultural',
        'diesel, burned in agricultural machinery'
        
    ]
    of_territory_input_new = [
        'nutrient supply from calcium nitrate'
    ]

    test_in_te_new = []
    test_of_te_new = []
    test_mx_te_new = []


    proc_name_value = process_name['name']

    # three energy crops have different input name
    # direct emissions
    for bio in process_name.biosphere():
        inpt_name_value = 'na'
        inpt_type_value = 'direct emission'
        inpt_amnt_value = 'na'

        emis_name_value = bio.as_dict()['name']
        emis_cate_value = bio.input.as_dict()['categories']
        emis_amnt_value = bio.as_dict()['amount']
        emis_sorc_value = 'in-territory emission'

        proc_name.append(proc_name_value)
        inpt_name.append(inpt_name_value)
        inpt_type.append(inpt_type_value)
        inpt_amnt.append(inpt_amnt_value)

        emis_name.append(emis_name_value)
        emis_sorc.append(emis_sorc_value)
        emis_cate.append(emis_cate_value)
        emis_amnt.append(emis_amnt_value)

        # technosphere
    for tech in process_name.technosphere():

        inpt_name_value = tech.input.as_dict()['name']

        inpt_amnt_value = tech['amount']

        # the inputs type in-off-mix information
        if tech['name'] in in_territory_input_new:
            test_in_te_new.append(tech['name'])
            inpt_type_value = 'in-territory input'

            for m in iw_method[:1]:

                lcia_name_value = m[2]

                lca = bw.LCA({tech.input: inpt_amnt_value}, m)
                lca.lci()
                lca.lcia()
                lcia_dfsc_value = lca.score

                for e in lca.biosphere_dict:
                    emis = bw.get_activity(e)
                    emis_name_value = emis.as_dict()['name']
                    emis_cate_value = emis.as_dict()['categories']
                    emis_sorc_value = 'in-territory emission'

                    emis_row = lca.biosphere_dict[emis]
                    emis_amnt_value = lca.inventory[emis_row, :].sum()

                    proc_name.append(proc_name_value)
                    inpt_name.append(inpt_name_value)
                    inpt_type.append(inpt_type_value)
                    inpt_amnt.append(inpt_amnt_value)

                    emis_name.append(emis_name_value)
                    emis_sorc.append(emis_sorc_value)
                    emis_cate.append(emis_cate_value)
                    emis_amnt.append(emis_amnt_value)

        elif tech['name'] in of_territory_input_new:
            test_of_te_new.append(tech['name'])
            inpt_type_value = 'off-territory input'
            for m in iw_method[:1]:

                lcia_name_value = m[2]

                lca = bw.LCA({tech.input: inpt_amnt_value}, m)
                lca.lci()
                lca.lcia()
                lcia_dfsc_value = lca.score

                for e in lca.biosphere_dict:
                    emis = bw.get_activity(e)
                    emis_name_value = emis.as_dict()['name']
                    emis_cate_value = emis.as_dict()['categories']
                    emis_sorc_value = 'off-territory emission'

                    emis_row = lca.biosphere_dict[emis]
                    emis_amnt_value = lca.inventory[emis_row, :].sum()

                    proc_name.append(proc_name_value)
                    inpt_name.append(inpt_name_value)
                    inpt_type.append(inpt_type_value)
                    inpt_amnt.append(inpt_amnt_value)

                    emis_name.append(emis_name_value)
                    emis_sorc.append(emis_sorc_value)
                    emis_cate.append(emis_cate_value)
                    emis_amnt.append(emis_amnt_value)

        else:
            test_mx_te_new.append(tech['name'])
            inpt_type_value = 'mix-territory input'

            for tech2 in tech.input.technosphere():

                tech_name_value = tech2.input.as_dict()['name']
                tech2_amnt = tech2['amount']

                for tech3 in tech2.input.technosphere():

                    tech3_amnt = tech3['amount']

                    tect_amnt = tech2_amnt * tech3_amnt

                    # calculate the LCA score for these technosphere input of the mix-input
                    for m in iw_method[:1]:
                        lcia_name_value = m[2]

                        lca = bw.LCA({tech3.input: inpt_amnt_value * tect_amnt}, m)
                        lca.lci()
                        lca.lcia()
                        lcia_dfsc_value = lca.score

                        for e in lca.biosphere_dict:
                            emis = bw.get_activity(e)
                            emis_name_value = emis.as_dict()['name']
                            emis_cate_value = emis.as_dict()['categories']
                            emis_sorc_value = 'off-territory emission'

                            emis_row = lca.biosphere_dict[emis]
                            emis_amnt_value = lca.inventory[emis_row, :].sum()

                            proc_name.append(proc_name_value)
                            inpt_name.append(inpt_name_value)
                            inpt_type.append(inpt_type_value)
                            inpt_amnt.append(inpt_amnt_value)

                            emis_name.append(emis_name_value)
                            emis_sorc.append(emis_sorc_value)
                            emis_cate.append(emis_cate_value)
                            emis_amnt.append(emis_amnt_value)

                for bio3 in tech2.input.biosphere():
                    emis_name_value = bio3.as_dict()['name']
                    emis_sorc_value = 'in-territory emission'
                    emis_cate_value = bio3.input.as_dict()['categories']
                    emis_amnt_value = inpt_amnt_value * tech2_amnt * bio3['amount']

                    proc_name.append(proc_name_value)
                    inpt_name.append(inpt_name_value)
                    inpt_type.append(inpt_type_value)
                    inpt_amnt.append(inpt_amnt_value)

                    emis_name.append(emis_name_value)
                    emis_sorc.append(emis_sorc_value)
                    emis_cate.append(emis_cate_value)
                    emis_amnt.append(emis_amnt_value)

    tef_new = pd.DataFrame()
    tef_new['proc_name'] = proc_name
    tef_new['inpt_name'] = inpt_name
    tef_new['inpt_type'] = inpt_type
    tef_new['inpt_amnt'] = inpt_amnt
    tef_new['emis_name'] = emis_name
    tef_new['emis_sorc'] = emis_sorc
    tef_new['emis_cate'] = emis_cate
    tef_new['emis_amnt'] = emis_amnt  # emission amount for each input/ha

    return tef_new


def overlay4(gdf1, gdf2, gdf3, gdf4):
    
    o1 = gpd.overlay(gdf1, gdf2, how='intersection')
    o2 = gpd.overlay(o1, gdf3, how='intersection')
    o3 = gpd.overlay(o2, gdf4, how='intersection')

    # o_sort = o3.sort_values(by=['cc_ID', 'fw_ID', 'yd_ID', 'Communes'])
    gdf_res = o3.reset_index().drop(columns='index')
    
    return gdf_res

def overlay7(gdf1, gdf2, gdf3, gdf4, gdf5, gdf6, gdf7):
    
    o1 = gpd.overlay(gdf1, gdf2, how='intersection')
    o2 = gpd.overlay(o1, gdf3, how='intersection')
    o3 = gpd.overlay(o2, gdf4, how='intersection')
    o4 = gpd.overlay(o3, gdf5, how='intersection')
    o5 = gpd.overlay(o4, gdf6, how='intersection')
    o6 = gpd.overlay(o5, gdf7, how='intersection')

    # o_sort = o3.sort_values(by=['cc_ID', 'fw_ID', 'yd_ID', 'Communes'])
    gdf_res = o6.reset_index().drop(columns='index')
    
    return gdf_res

def cal_ef_truck(inventory, name, iw_method):

    inpt_name = []
    inpt_type = []  # mix-territory input
    inpt_amnt = []

    emis_name = []
    emis_sorc = []  # define if it is in-territory emission or off-territory emission
    #                 'in-territory emission', 'off-territory emission'
    emis_cate = []
    emis_amnt = []

    # direct emissions
    for bio in inventory.biosphere():
        inpt_name_value = 'na'
        inpt_type_value = 'direct emission'
        inpt_amnt_value = 'na'

        emis_name_value = bio.as_dict()['name']
        emis_cate_value = bio.input.as_dict()['categories']
        emis_amnt_value = bio.as_dict()['amount']  # calculate the amount for 1ton kilometer of product
        emis_sorc_value = 'in-territory emission'

        inpt_name.append(inpt_name_value)
        inpt_type.append(inpt_type_value)
        inpt_amnt.append(inpt_amnt_value)

        emis_name.append(emis_name_value)
        emis_sorc.append(emis_sorc_value)
        emis_cate.append(emis_cate_value)
        emis_amnt.append(emis_amnt_value)

        # technosphere
    for tech in inventory.technosphere():

        inpt_name_value = tech.input.as_dict()['name']

        inpt_type_value = 'off-territory input'
        inpt_amnt_value = tech['amount']

        print(inpt_type_value, inpt_name_value)


        # calculate the LCA score for these technosphere input of the off-input
        for m in iw_method[:1]:
            lcia_name_value = m[2]

            lca = bw.LCA({tech.input: inpt_amnt_value}, m)
            lca.lci()
            lca.lcia()
            lcia_dfsc_value = lca.score

            for e in lca.biosphere_dict:
                emis = bw.get_activity(e)
                emis_name_value = emis.as_dict()['name']
                emis_cate_value = emis.as_dict()['categories']
                emis_sorc_value = 'off-territory emission'

                emis_row = lca.biosphere_dict[emis]
                emis_amnt_value = lca.inventory[emis_row, :].sum()

                inpt_name.append(inpt_name_value)
                inpt_type.append(inpt_type_value)
                inpt_amnt.append(inpt_amnt_value)

                emis_name.append(emis_name_value)
                emis_sorc.append(emis_sorc_value)
                emis_cate.append(emis_cate_value)
                emis_amnt.append(emis_amnt_value)


    emission = pd.DataFrame()
    emission['inve_name'] = name
    emission['inpt_name'] = inpt_name
    emission['inpt_amnt'] = inpt_amnt
    emission['emis_name'] = emis_name
    emission['emis_sorc'] = emis_sorc
    emission['emis_cate'] = emis_cate
    emission['emis_amnt'] = emis_amnt  # emission amount for each input/ton kilometer

    return emission

    