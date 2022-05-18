# coding: utf-8

### Filename: driver.py
### Language: Python v2.7
### Created by: Jesse Gourevitch

import os
import gc
import ogr
import itertools
import numpy as np
import pandas as pd

from time import time, ctime
from scipy import interpolate
from datetime import datetime

import constants, paths, preprocessing, utils

### Ignore all numpy errors
np.seterr(all='ignore')

def execute(ws, huc12_list, n, resolutions, percentiles,reach_type):
    """Execute model functions, depending on boolean operators in driver_args.

    Parameters:

    Returns:
        None
    """

    ### Generate RVS for parameters
    rvs_dict = utils.generate_rvs(n)

    ### Identify missing reaches
    print_missingreaches = False
    if print_missingreaches == True:
        utils.streamstats_qaqc()

    ### Get watershed ID from dictionary
    ws_id = constants.ws_id_dict[ws]

    ### Iterate through HUC12s
    for huc12 in huc12_list:
        huc12_id = huc12.split('_')[1]
        t0 = datetime.now()
        t0_str = t0.strftime('%H:%M:%S')

        ### Initialize empty Pandas DataFrames
        rc_df, RIstage_df = [pd.DataFrame()]*2

        print('')
        print('*********************************')
        print('Processing %s %s @ %s' %(ws_id, huc12_id, t0_str))
        print('*********************************')

        ### Iterate through spatial resolutions
        for cellsize in resolutions:

            ### Get paths_dict from paths.py
            ### ws ==> HUC8 long-form name (IE Winooski_River)
            ### ws_id ==> HUC8 abbreviation (IE WIN)
            paths_dict = paths.gen_paths(ws, ws_id, huc12_id,reach_type)


            if constants.only_generate_HAND:
                
                preprocess_data(paths_dict, cellsize)
            else:
                ### Preprocess data
                paths_dict = preprocess_data(paths_dict, cellsize)

                ### Read data files to memory
                data_dict = preprocessing.read_data2mem(paths_dict, cellsize)
                
                ### Execute Monte Carlo Simulation
                rc_df, RIstage_df = monte_carlo(
                    paths_dict, data_dict, rvs_dict, n, huc12_id, ws_id,
                    cellsize, rc_df, RIstage_df)
                
                export_logbooks(rc_df, RIstage_df, paths_dict)
    
                for p in percentiles:        
                    map_inun(RIstage_df, data_dict, p, cellsize, paths_dict)
    
                ### Empty garbage can
                data_dict = None
                gc.collect()

            ### Print elapsed time
            utils.print_elapsedtime(t0)
                
    return None


def preprocess_data(paths_dict, cellsize):

    ### Generate HAND layer
    breakhere=1
    preprocessing.generate_hand(paths_dict)
    
    ### Resample HAND and LULC layers to a lower resolution
    paths_dict = preprocessing.resample_data(paths_dict, cellsize)

    ### For 1m resolution HAND layers, impute missing values 
    paths_dict = preprocessing.impute_hand(paths_dict, cellsize)

    ### Generate Theissen raster
    preprocessing.generate_theissen(paths_dict, cellsize)

    ### Extract DEM values to stream points
    preprocessing.extract_dem2points(paths_dict)

    return paths_dict


def monte_carlo(paths_dict, data_dict, rvs_dict, n, huc12_id, 
                ws_id, cellsize, rc_df, RIstage_df):

    ### Initialize reach list from stream stats
    ss_df = data_dict['streamstats_df']
    ss_df2 = ss_df[ss_df['Huc12']=='%s_%s' %(ws_id, huc12_id)]
    reach_list = pd.Series(ss_df2['Name'].unique())

    stream_df = data_dict['stream_df']
    
    #Check that the reaches consist of a single polyline
    max_segments_per_reach = stream_df.groupby(['Code']).count()['Length'].max()
    if max_segments_per_reach != 1:
        print("******\n   ERROR: a reach segment consists of {} polylines - exiting\n******".format(max_segments_per_reach))
        return 
    
    reach_list = list(reach_list[reach_list.isin(stream_df['Code'])])       

    ### Iterate through simulations
    for i in range(n):
        params_dict = utils.rvs2dict(i, rvs_dict)

        if i%(n/10) == 0: 
            print ('\t%s - Starting iteration #%d' %(ctime(), i))

        if i==0:
            rc_df = generate_base_rc(i, rc_df,   data_dict, reach_list, cellsize)
        else:
            rc_df = generate_rc_frombase(i, rc_df, params_dict, reach_list, cellsize)

        ### Get stage for each reach for each recurrance interval
        for yr, reach in itertools.product(constants.interval_list, reach_list):
            reach_stage_dict = get_reach_stage(
                i, huc12_id, ws_id, rc_df, data_dict, 
                params_dict, yr, cellsize, reach)
            RIstage_df = RIstage_df.append(reach_stage_dict, ignore_index=True)
    return rc_df, RIstage_df


def generate_base_rc(i, rc_df, data_dict, reach_list, cellsize):

    ### Create a copy of the lulc array
    lulc_array = data_dict['lulc_array'].copy()

    ### Replace LULC codes with Manning's n values
    lulc2mannings_dict = constants.lulc2mannings_dict
    n_array = np.zeros(lulc_array.shape)
    for lulc_code in lulc2mannings_dict.keys():
        n_array[lulc_array==int(lulc_code)] = (
            lulc2mannings_dict[lulc_code])

    ### Initialize data from data dict
    hand_array = data_dict['hand_array']
    thiessen_array = data_dict['thiessen_array']
    streampoints_df = data_dict['streampoints_df']
    stream_df = data_dict['stream_df']

    ### Iterate over through user-specified stages and over reaches in HUC12
    for stage, reach in itertools.product(constants.stage_list, reach_list):
        reach = float(reach)

        ### Subset HAND array to where HAND is <= to stage 
        mask = (hand_array<=stage) & (thiessen_array==reach)
        hand_array_masked = hand_array[mask]

        ### Calculate volume & surface area
        volume = np.nansum(stage - hand_array_masked) * (cellsize**2)
        surface_area = (np.count_nonzero(~np.isnan(hand_array_masked)) * 
                       (cellsize**2))
        
        ### Calculate cross-sectional area
        length = float(stream_df['Length'][stream_df['Code']==reach])
        # length = float(stream_df['SWLENGTH'][stream_df['Code']==reach])
        xs_area = volume / length
        
        ### Calculate hydraulic radius
        hydraulic_radius = volume / surface_area
        
        ### Calculate average Manning's coefficient, weighted by volume
        n_array_masked = n_array[mask]
        n_volume_total = 0

        for n in np.unique(n_array_masked[n_array_masked!=0]):
            hand_subset = hand_array_masked[n_array_masked==n]

            n_volume = np.nansum(stage - hand_subset) * cellsize**2
            n_volume_total += (n_volume * n)

        mannings_mean = n_volume_total / volume

        if mannings_mean == 0:
            mannings_mean = 0.089666666666 # Avg of all Manning's values
        
        ### Calculate slope based on DEM min and max
        streampoints_df2 = streampoints_df[
            streampoints_df['Code']==reach]
        elev_min = streampoints_df2['ELEV'].min()
        elev_max = streampoints_df2['ELEV'].max()
        slope = (elev_max - elev_min) / length

        ### If calculated slope is < 0.00001, use NHD+ slope instead
        if slope < 0.00001:
            slope = float(stream_df['SloNHDcalc']
                [stream_df['Code']==reach])

            ### If NHD+ slope is still less than 0.00001, set to 0.00001
            if slope < 0.00001: slope = 0.00001

        ### Calculate Q using Mannings equation
        Q = ((1 / mannings_mean) * np.power(hydraulic_radius, (2.0/3.0)) * 
             (np.sqrt(slope)) * (xs_area))
        if volume == 0: Q = 0

        rc_row_dict = {
            'N_SIM': int(i),
            'RESOLUTION': cellsize, 
            'REACH': int(reach),
            'STAGE': stage,
            'LENGTH': length,
            'VOLUME': volume,
            'SA': surface_area,
            'XS_AREA': xs_area,
            'H_RADIUS': hydraulic_radius,
            'MANNINGS': mannings_mean,
            'SLOPE': slope,
            'DISCHARGE': Q
            }
        rc_df = rc_df.append(rc_row_dict, ignore_index=True)

    return rc_df


def generate_rc_frombase(i, rc_df, params_dict, reach_list, cellsize):
    ### Iterate over through user-specified stages and over reaches in HUC12
    for stage, reach in itertools.product(constants.stage_list, reach_list):
        reach = float(reach)

        ### Get baseline reach parameters
        rc_df_subset = rc_df[(rc_df['RESOLUTION']==cellsize) &
                             (rc_df['N_SIM']==0) &
                             (rc_df['REACH']==reach) &
                             (rc_df['STAGE']==stage)]

        xs_area = float(rc_df_subset['XS_AREA'])
        hydraulic_radius = float(rc_df_subset['H_RADIUS'])
        mannings_mean = float(rc_df_subset['MANNINGS'])
        slope = float(rc_df_subset['SLOPE'])
        volume = float(rc_df_subset['VOLUME'])
        length = float(rc_df_subset['LENGTH'])
        surface_area = float(rc_df_subset['SA'])

        ### Multiple calculated reach parameters by RVS scalars
        mannings_mean *= params_dict['mannings_rvs']
        if slope < 0.001:
            xs_area *= params_dict['low_grade_xsarea_rvs']
            hydraulic_radius *= params_dict['low_grade_xsarea_rvs']
            slope *= params_dict['low_grade_slope_rvs']
        else:
            xs_area *= params_dict['high_grade_xsarea_rvs']
            hydraulic_radius *= params_dict['high_grade_xsarea_rvs']
            slope *= params_dict['high_grade_slope_rvs']

        ### If slope is less than 0.00001, set to 0.00001
        if slope < 0.00001: slope = 0.00001

        ### Calculate Q using Mannings equation
        Q = ((1 / mannings_mean) * np.power(hydraulic_radius, (2.0/3.0)) * 
             (np.sqrt(slope)) * (xs_area))
        if volume == 0: Q = 0

        rc_row_dict = {
            'N_SIM': int(i),
            'RESOLUTION': cellsize, 
            'REACH': int(reach),
            'STAGE': stage,
            'LENGTH': length,
            'VOLUME': volume,
            'SA': surface_area,
            'XS_AREA': xs_area,
            'H_RADIUS': hydraulic_radius,
            'MANNINGS': mannings_mean,
            'SLOPE': slope,
            'DISCHARGE': Q
            }
        rc_df = rc_df.append(rc_row_dict, ignore_index=True)
    
    return rc_df


def get_reach_stage(i, huc12_id, ws_id, rc_df, data_dict, 
                    params_dict, yr, cellsize, reach):
    reach = float(reach)
    ss_df = data_dict['streamstats_df']
    Q = ss_df['Value'][(ss_df['Huc12']=='%s_%s' %(ws_id, huc12_id)) & 
                         (ss_df['Name']==reach) &
                         (ss_df['StatName']=='%s Year Peak Flood' %str(yr))]
    if len(Q) == 0: Q = 0
    else: Q = float(Q)

    ### Convert Q from (ft^3 / sec) to (m^3 / sec)
    Q *= 0.0283168

    ### Apply uncertainty to Q
    Q *= params_dict['discharge_rvs_Q%d' %yr]

    rc_df2 = rc_df[(rc_df['RESOLUTION']==cellsize) & 
                   (rc_df['N_SIM']==int(i)) &
                   (rc_df['REACH']==reach)]

    ### Linear interpolation between discharge values
    Q_low = (rc_df2['DISCHARGE'][rc_df2['DISCHARGE'] < Q]).max()
    Q_high = (rc_df2['DISCHARGE'][rc_df2['DISCHARGE'] > Q]).min()

    if Q_low > 0:
        stage_low = rc_df2['STAGE'][rc_df2['DISCHARGE'] == Q_low]
        stage_low = float(stage_low)
    else:
        stage_low = np.nan
    
    if Q_high > 0:
        stage_high = rc_df2['STAGE'][rc_df2['DISCHARGE'] == Q_high]
        stage_high = float(stage_high)

    else: 
        stage_high = np.nan

    if np.isnan(stage_low) or np.isnan(stage_high):
        stage = 0
    else:
        f = interpolate.interp1d([Q_low, Q_high], [stage_low, stage_high])
        stage = f(Q)

    reach_stage_dict ={
        'REACH': reach,
        'STAGE': stage,
        'Q': Q,
        'RI': yr,
        'N_SIM': i,
        'RESOLUTION': cellsize
        }

    return reach_stage_dict


def map_inun(RIstage_df, data_dict, p, cellsize, paths_dict):

    print ('\t%s - Mapping inundation for %dth percentile' %(ctime(), p))

    ### Initialize data from data_dict
    hand_array = data_dict['hand_array']
    hand_gt = data_dict['hand_gt']
    thiessen_array = data_dict['thiessen_array']
    inun_array = np.where(~np.isnan(hand_array), 0, np.nan)

    inun_shp_list = []
    stage_list = []
    for yr in constants.interval_list:
        for reach in RIstage_df['REACH'].unique():
            reach = float(reach)
            stage = RIstage_df['STAGE'][(RIstage_df['RI']==yr) & 
                                        (RIstage_df['REACH']==reach) & 
                                        (RIstage_df['RESOLUTION']==cellsize)]
            stage_p = np.percentile(stage, 100-p)
            stage_list.append(stage_p)

            ### Fill in values for inundation array
            mask = (hand_array <= stage_p) & (thiessen_array==reach)
            depth_array = stage_p - hand_array[mask]
            inun_array[mask] = depth_array
                
            ### Convert zero values to NaN
            inun_array = np.where(inun_array==0, np.nan, inun_array)

        ### Export rasters
        inun_raster_uri = export_raster(inun_array, paths_dict, yr, p, cellsize)
        
    return None


def export_raster(inun_array, paths_dict, yr, p, cellsize):
    inun_rasters_dir = paths_dict['inundation_rasters_dir']
    
    hand_uri = paths_dict['hand_uri']

    if np.all(np.isnan(inun_array)) == True:
        print ('\tERROR: All px in Q%d inun_array are NaN' %yr)
        print ('\tERROR (cont): Skipping inun_array export')
        print('Try checking that the "HUC12", "Name", and "StatName" fields')
        print('\tare populated correctly in the streamstats csv file')
        quit()
    else:
        inun_raster_uri = os.path.join(
            inun_rasters_dir, 'Q%d_Depth_%dm_p%d.tif' %(yr, cellsize, p))
        utils.array2raster(inun_array, hand_uri, inun_raster_uri)

    return inun_raster_uri


def export_logbooks(rc_df, RIstage_df, paths_dict):
    print ('\t%s - Exporting logbooks' %ctime())

    ### Set order of columns in rc_df and RIstage_df
    rc_df = rc_df[[
        'RESOLUTION', 'N_SIM', 'REACH', 'STAGE', 
        'LENGTH', 'SA', 'VOLUME', 'XS_AREA', 'H_RADIUS', 
        'MANNINGS', 'SLOPE', 'DISCHARGE']]

    RIstage_df = RIstage_df[['RESOLUTION', 'N_SIM', 
                             'REACH', 'RI', 'Q', 'STAGE']]

    ### Sort values in rc_df and RIstage_df
    rc_df = rc_df.sort_values(
        by=['RESOLUTION', 'N_SIM', 'REACH', 'STAGE'])
    RIstage_df = RIstage_df.sort_values(by=['N_SIM', 'REACH', 'RI'])

    ### Export logbooks to CSVs
    rc_logbook_fn = 'rc_logbook.csv'
    rc_logbook_uri = os.path.join(
        paths_dict['logbooks_dir'], rc_logbook_fn)
    rc_df.to_csv(rc_logbook_uri, index=False)

    RIstage_logbook_fn = 'RIstage_logbook.csv'
    RIstage_logbook_uri = os.path.join(
        paths_dict['logbooks_dir'], RIstage_logbook_fn)
    RIstage_df.to_csv(RIstage_logbook_uri, index=False)

    return None

