# coding: utf-8

### Filename: paths.py
### Language: Python v2.7
### Created by: Jesse Gourevitch

import os
import platform

import constants


def gen_paths(ws, ws_id, huc12_id, reach_type, data_folder):
    """Generate paths to directories, sub-directories, data input and 
        output files, and parameters.
    
    Parameters:
        ws (str): watershed of interest (e.g. 'Otter_Creek')

    Returns:
        paths_dict (dict): dictionary containing data paths and parameter values
    """
    
    
    '''
    ws_id = 'WIN'
    huc12_id = '0501'
    ws = 'Winooski_River'
    reach_type = 'SGA'
    stem = 'D:\\Github\\mad_river'
    '''

    ### Initialize path to data directory
    if not data_folder:
        stem = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    else:
        stem = data_folder
    data_dir = os.path.join(stem, 'Data')

    huc8_code = constants.huc8_code_dict[ws_id]
    huc12_code = huc8_code + huc12_id
    
    threshold_flow = int(round(constants.threshold_flow,0))
    
    cd = os.getcwd()
    if cd.split(os.path.sep)[-1] == 'source':
        alg_type = 'static'
    else:
        alg_type = 'dynamic'

    ### Create paths dictionary
    paths_dict = {
        ### Initialize path to TauDEM toolbox
        'taudem_toolbox_path': 'C:/Program Files/TauDEM/TauDEM5Arc/TauDEM Tools.tbx',
        
        ### Store HUC information
        'ws_id':ws_id,
        'huc12_id':huc12_id,
        'reach_type':reach_type,
        'threshold':threshold_flow,

        ### Initialize paths to sub-directories
        'dem_dir':os.path.join(data_dir, ws,'Buffered_DEM'),
        'dem_derivates_dir': os.path.join(data_dir, ws, 'DEM_Derivatives', f'thresh_{threshold_flow}', f'HUC12_{huc12_id}'),
        'streampoints_dir': os.path.join(data_dir, ws, 'Stream_Points', reach_type, f'HUC12_{huc12_id}'),
        'thiessen_dir': os.path.join(data_dir, ws, 'Thiessen_Polygons', reach_type, f'HUC12_{huc12_id}'),
        'inundation_rasters_dir': os.path.join(data_dir, ws, 'Inundation_Rasters', alg_type, reach_type, f'thresh_{threshold_flow}', f'HUC12_{huc12_id}'),
        'inun_shp_dir': os.path.join(data_dir, ws, 'Inundation_Shapefiles', alg_type, reach_type, f'thresh_{threshold_flow}', f'HUC12_{huc12_id}'),
        'logbooks_dir': os.path.join(data_dir, ws, 'Output_Logbooks', alg_type, reach_type, f'thresh_{threshold_flow}', f'HUC12_{huc12_id}'),
        #path to the DEM covering the entire HUC8
        'huc8_DEM_dir': os.path.join(data_dir, ws, 'Buffered_DEM', 'huc8_DEM'),
        #path to the huc8-level stream raster
        # 'huc8_stream_dir': os.path.join(data_dir, ws, 'DEM_Derivatives', 'HUC8_derivatives'),
        #Location to store the intermediate steps in generating the huc8 scale
        #flow accumulation raster
        'huc8_derivatives_dir': os.path.join(data_dir, ws, 'DEM_Derivatives', 'HUC8_derivatives'),
        'temp_flow_acc_dir': os.path.join(data_dir, ws, 'DEM_Derivatives', 'HUC8_derivatives'),
        'huc8_streamline_dir': os.path.join(data_dir, ws, 'DEM_Derivatives', 'HUC8_derivatives'),
        #huc8-level pit-filled DEM
        'huc8_dem_buff_filled_uri': os.path.join(data_dir, ws, 'DEM_Derivatives', 'HUC8_derivatives', 'temporary_flow_accumulation_files', f'{ws_id}_pitfilled.tif'),
        
        #HAND generation logs
        'HAND_log_uri': os.path.join(data_dir, ws, 'HAND_generation_logs', f'{ws_id}_{huc12_id}.txt'),
        'HAND_log_dir': os.path.join(data_dir, ws, 'HAND_generation_logs'),
        
        #headwater mask directory
        # 'headwater_mask_dir': os.path.join(data_dir, ws, 'DEM_Derivatives', 'HUC8_derivatives','headwater_mask'),
        'watershed_output_uri': os.path.join(data_dir,ws,'Watershed','watershed_output.shp'),
        
        
        
        'hand_dir': os.path.join(data_dir, ws, 'HAND'),
        
        
        'hash_dir': os.path.join(data_dir, ws, 'File_Hashes'),

        ### Initialize paths to data inputs & outputs
        'dem_buff_uri': os.path.join(data_dir, ws, 'Buffered_DEM', f'{ws_id}_{huc12_id}.tif'),
        'ws_buff_uri': os.path.join(data_dir, ws, 'Watershed_Buffer', f'{ws_id}_{huc12_id}_BF.shp'),
        'ws_uri': os.path.join(data_dir, ws, 'Watershed', f'{ws_id}_{huc12_id}.shp'),
        'lulc_uri': os.path.join(data_dir, ws, 'LULC', f'{ws_id}_{huc12_id}_LC.tif'),
        'hand_huc8_5m_uri': os.path.join(data_dir, f'LCB/HAND/{ws_id}_HAND_5m.tif'),
        'hand_lcb_7point4m_uri': os.path.join(data_dir, 'LCB/HAND/LCB_HAND_7point4m.tif'),
        #path to the DEM covering the entire HUC8
        'huc8_dem_uri': os.path.join(data_dir, ws, 'Buffered_DEM', 'huc8_DEM', f'{ws}.tif'),
        'hand_uri': os.path.join(data_dir, ws, 'HAND', f'{ws_id}_{huc12_id}_thresh-{threshold_flow}_HAND.tif'),
        #huc8-level stream raster from flow accumulation calculation
        'stream_huc8_uri': os.path.join(data_dir, ws, 'DEM_Derivatives', 'HUC8_derivatives', f'{ws_id}_thresh-{threshold_flow}_stream_network.tif'),
        
        #reach type dependant uri variables
        'stream_uri': os.path.join(data_dir, ws, 'Stream_Polylines', reach_type, f'{ws_id}_{huc12_id}_{reach_type}.shp'),
        'ratingcurve_uri': os.path.join(data_dir, ws, 'Rating_Curves', alg_type, reach_type, f'thresh_{threshold_flow}', f'{ws_id}_{huc12_id}_{reach_type}_RC.csv'),
        'streamstats_uri': os.path.join(data_dir, ws, 'Stream_Stats', f'Stream_Stats_{reach_type}.csv')
        }

    ### If directory does not exist, create it
    paths = list(paths_dict.keys())[1:]
    
    for path in paths:
        
        if (not os.path.exists(paths_dict[path])) and ('dir' in path):
            breakhere = 1
            os.makedirs(paths_dict[path])

    return paths_dict



