import numpy as np
import pandas as pd
import geopandas as gpd
import os
import matplotlib.pyplot as plt

BASIN_DIRECTORY = r'D:\HAND\Data\Winooski_River'



def condense_logbooks(logbook_dir, hydrology_dir, export_path):
    # Get list of paths to mine
    path_list = list()
    for folder in os.listdir(logbook_dir):
        tmp_folder = os.path.join(logbook_dir, folder)
        if 'rc_logbook.csv' in os.listdir(tmp_folder) and 'RIstage_logbook.csv' in os.listdir(tmp_folder):
            path_list.append(
                (os.path.join(tmp_folder, 'rc_logbook.csv'), os.path.join(tmp_folder, 'RIstage_logbook.csv')))
    print(f'Found {len(path_list)} basins to condense logbooks for')

    # Condense logbooks
    working_rc = pd.DataFrame()
    working_ri = pd.DataFrame()
    counter = 1
    for tmp_paths in path_list:
        print(f' - {counter} / {len(path_list)}', end='\r')
        counter += 1

        rc_logbook = pd.read_csv(tmp_paths[0])
        ri_logbook = pd.read_csv(tmp_paths[1])
        base_rc = rc_logbook[rc_logbook['N_SIM'] == 0]
        base_ri = ri_logbook[ri_logbook['N_SIM'] == 0]

        if len(working_rc) == 0:
            working_rc = base_rc
            working_ri = base_ri
        else:
            working_rc = pd.concat([working_rc, base_rc])
            working_ri = pd.concat([working_ri, base_ri])
    print('\n')

    # Append drainage area information
    hydrology_df = pd.read_csv(hydrology_dir)
    hydrology_df = hydrology_df[['Name', 'AreaSqMi']].drop_duplicates()
    working_ri = pd.merge(working_ri, hydrology_df, how='left', left_on='REACH', right_on='Name')

    # Export
    export_paths = (os.path.join(export_path, 'combo_rc_logbook.csv'), os.path.join(export_path, 'combo_ri_logbook.csv'))
    working_rc.to_csv(export_paths[0])
    working_ri.to_csv(export_paths[1])
    print(f'exported {export_paths[0]} and {export_paths[1]}')
    return export_paths


def export_all_curves(rc_path, ri_path, export_dir, thresholds=None):
    # Load data
    rc_df = pd.read_csv(rc_path)
    ri_df = pd.read_csv(ri_path)

    # Extract reach list and set up processed data dictionaries
    reaches = rc_df['REACH'].unique()
    print(f'Found {len(reaches)} reaches to process')
    intervals = ri_df['RI'].unique()
    raw_out_dict = {'REACH': list(), 'QNORM': list(), 'CH_SSP': list(), 'OB_SSP': list()}
    interpolated_out_dict = {'REACH': list()}
    for ri in intervals:
        interpolated_out_dict[f'CH_SSP_Q{ri}'] = list()
        interpolated_out_dict[f'OB_SSP_Q{ri}'] = list()

    counter = 1
    for reach in reaches:
        print(f' - {counter} / {len(reaches)}', end='\r')
        counter += 1

        reach_filter_rc = rc_df['REACH'] == reach
        reach_filter_ri = ri_df['REACH'] == reach

        q2_flow = ri_df[reach_filter_ri & (ri_df['RI'] == 2)]['Q'].item()
        drainage_area = float(ri_df[reach_filter_ri & (ri_df['RI'] == 2)]['AreaSqMi'].item())

        # Import base RC
        flows = np.array(rc_df[reach_filter_rc]['DISCHARGE'])
        stages = np.array(rc_df[reach_filter_rc]['STAGE'])
        sa = np.array(rc_df[reach_filter_rc]['SA'])
        volumes = np.array(rc_df[reach_filter_rc]['VOLUME'])
        mannings = np.array(rc_df[reach_filter_rc]['MANNINGS'])
        length = np.array(rc_df[reach_filter_rc]['LENGTH'].unique()[0])
        slope = np.array(rc_df[reach_filter_rc]['SLOPE'].unique()[0])

        # Split channel and overbank
        threshold = 0
        for t in thresholds:
            if (drainage_area >= thresholds[t][0]) & (drainage_area < thresholds[t][1]):
                threshold = t

        sa_at_split = np.interp(threshold, stages, sa)
        volume_at_split = np.interp(threshold, stages, volumes)

        ch_sa = np.array([a if a <= sa_at_split else sa_at_split for a in sa])
        ch_v = np.array([v if v <= volume_at_split else volume_at_split + (sa_at_split * (s - threshold)) for v, s in zip(volumes, stages)])
        ob_sa = sa - ch_sa
        ob_v = volumes - ch_v

        # Hydraulic modeling
        ch_q = ((ch_v / length) * ((ch_v / ch_sa) ** (2 / 3)) * (slope ** 0.5)) / mannings
        ob_q = ((ob_v / length) * ((ob_v / ob_sa) ** (2 / 3)) * (slope ** 0.5)) / mannings
        np.nan_to_num(ob_q, copy=False)
        combo_q = ch_q + ob_q
        qnorm = combo_q / q2_flow

        ch_ssp = (9810 * ch_q * slope) / (ch_sa / length)
        ob_ssp = (9810 * ob_q * slope) / (ob_sa / length)
        np.nan_to_num(ob_ssp, copy=False)

        # Interpolate SSP at RIs
        interpolated_out_dict['REACH'].append(reach)
        for ri in intervals:
            tmp_flow = ri_df[reach_filter_ri & (ri_df['RI'] == ri)]['Q'].item()
            interpolated_out_dict[f'CH_SSP_Q{ri}'].append(np.interp(tmp_flow, combo_q, ch_q))
            interpolated_out_dict[f'OB_SSP_Q{ri}'].append(np.interp(tmp_flow, combo_q, ob_q))

        # Log data
        raw_out_dict['REACH'].extend([reach] * len(stages))
        raw_out_dict['QNORM'].extend(qnorm)
        raw_out_dict['CH_SSP'].extend(ch_ssp)
        raw_out_dict['OB_SSP'].extend(ob_ssp)
    print('\n')

    raw_out_df = pd.DataFrame(raw_out_dict)
    raw_out_path = os.path.join(export_dir, 'raw_ssp_curves.csv')
    raw_out_df.to_csv(raw_out_path)

    interpolated_out_df = pd.DataFrame(interpolated_out_dict)
    interpolated_out_path = os.path.join(export_dir, 'interpolated_ssp_curves.csv')
    interpolated_out_df.to_csv(interpolated_out_path)

    return raw_out_path, interpolated_out_path


def merge_with_gis(ssp_data_path, gis_dir, export_dir):
    # Merge all polylines shapefiles into one gdf
    paths = [os.path.join(gis_dir, f) for f in os.listdir(gis_dir) if f[-3:] == "shp"]
    print(f'Found {len(paths)} shapefiles to merge')
    working_gdf = gpd.GeoDataFrame(pd.concat([gpd.read_file(p) for p in paths],
                                             ignore_index=True), crs=gpd.read_file(paths[0]).crs)
    working_gdf = working_gdf[['Code', 'geometry']]
    working_gdf['Code'] = working_gdf['Code'].astype(int)

    # Join SSP data
    print('Joining SSP data')
    ssp_df = pd.read_csv(ssp_data_path)
    working_gdf = working_gdf.merge(ssp_df, how='left', left_on='Code', right_on='REACH')

    # Export
    export_path = os.path.join(export_dir, 'basin_ssp.shp')
    working_gdf.to_file(export_path)
    print(f'Exported to {export_path}')


def extract_ssp(basin_driectory, thresholds=None, reach_type='NHD', model_res='5179976'):
    ### Convenience method for user to run all SSP mining scripts ###

    # Set up ssp working directory
    working_dir = os.path.join(basin_driectory, 'Post-processing', 'ssp')
    if not os.path.exists(working_dir):
        if not os.path.exists(os.path.join(basin_driectory, 'Post-processing')):
            os.mkdir(os.path.join(basin_driectory, 'Post-processing'))
        os.mkdir(working_dir)

    # Combine logbooks for easier access and processing
    combo_rc_path, combo_ri_path = condense_logbooks(
        os.path.join(basin_driectory, 'Output_Logbooks', 'static', reach_type, f'thresh_{model_res}'),
        os.path.join(basin_driectory, 'Stream_Stats', f'Stream_Stats_{reach_type}.csv'),
        working_dir)

    # Extract SSP data
    raw_curves, interp_curves = export_all_curves(combo_rc_path, combo_ri_path, working_dir, thresholds=thresholds)

    # Merge SSP data with geospatial data (GIS polylines)
    merge_with_gis(interp_curves, os.path.join(basin_driectory, 'Stream_Polylines', 'NHD'), working_dir)



if __name__ == '__main__':
    threshold_dict = {0.5: [0, 2000],
                      1.0: [2000, 10000]}
    extract_ssp(BASIN_DIRECTORY, thresholds=threshold_dict)