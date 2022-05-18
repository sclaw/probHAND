import csv
import itertools
import os
import pandas as pd
import geopandas as gpd
from osgeo import gdal
import math
import numpy as np

TMP_HARDCODED_TREE_PATH = r"D:\Academic\UVM\ffi_floodplain_mapping_fall_2021\winooski_river\WIN_Tree.csv"
BASIN_IDS = ['WIN_0101', 'WIN_0102', 'WIN_0103', 'WIN_0201', 'WIN_0202', 'WIN_0203', 'WIN_0204', 'WIN_0301', 'WIN_0302', 'WIN_0401', 'WIN_0402', 'WIN_0403', 'WIN_0501', 'WIN_0502', 'WIN_0503', 'WIN_0504', 'WIN_0601', 'WIN_0602', 'WIN_0603', 'WIN_0604', 'WIN_0701', 'WIN_0702', 'WIN_0703', 'WIN_0704']
#BASIN_IDS = ["OTR_0101","OTR_0102","OTR_0103","OTR_0104","OTR_0105","OTR_0106","OTR_0107","OTR_0108","OTR_0109","OTR_0201","OTR_0202","OTR_0203","OTR_0301","OTR_0302","OTR_0303","OTR_0304","OTR_0305","OTR_0306","OTR_0307","OTR_0401","OTR_0402","OTR_0501","OTR_0502"]
#BASIN_IDS = ['MSQ_0101', 'MSQ_0102', 'MSQ_0103', 'MSQ_0104', 'MSQ_0105', 'MSQ_0201', 'MSQ_0202', 'MSQ_0203', 'MSQ_0204', 'MSQ_0301', 'MSQ_0302', 'MSQ_0401', 'MSQ_0402', 'MSQ_0501', 'MSQ_0502', 'MSQ_0503', 'MSQ_0601', 'MSQ_0602', 'MSQ_0603']
INTERVAL_LIST = [2, 5, 10, 25, 50, 100, 200, 500]
#PERCENTILES = [5, 10, 25, 50, 75, 90, 95]
PERCENTILES = [50]
RESOLUTION = 1
THRESHOLD = '5179976'
RIVER = 'Winooski_River'
REACH_TYPE = 'NHD'
DATA_PATH = r'D:\Academic\UVM\ffi_floodplain_mapping_fall_2021\winooski_river'


class Raster:

    def __init__(self, src_path):
        self.driver = gdal.GetDriverByName('GTiff')

        self.dataset = gdal.Open(src_path)
        self.band = self.dataset.GetRasterBand(1)
        self.band.ComputeStatistics(0)
        self.data = self.band.ReadAsArray()

        self.cols = self.dataset.RasterXSize
        self.rows = self.dataset.RasterYSize
        self.crs = self.dataset.GetProjectionRef()
        transform = self.dataset.GetGeoTransform()
        self.origin_x = transform[0]
        self.origin_y = transform[3]
        self.pixel_width = transform[1]
        self.pixel_height = transform[5]
        self.stats = self.band.GetStatistics(True, True)

    def index_to_geographic_pt(self, col, row):
        x = self.origin_x + (col*self.pixel_width) + (self.pixel_width/2)
        y = self.origin_y + (row*self.pixel_height) + (self.pixel_height/2)
        return x, y

    def geographic_pt_to_index(self, x, y):
        col = math.floor((x - self.origin_x) / self.pixel_width)
        row = math.floor((y - self.origin_y) / self.pixel_height)
        return col, row


def get_reach_tree(path):
    full_tree = {}  # relates parents to children
    parent_dict = {}
    with open(path, mode='r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            if len(row['Code']) != 0:
                full_tree[row['Code']] = {'parent': row['Parent'], 'basin': row['HUC12']}
                if row['Parent'] not in parent_dict:
                    parent_dict[row['Parent']] = [row['Code']]
                else:
                    parent_dict[row['Parent']].append(row['Code'])
    # remove all parents with one child.  remove all reaches not relevant to current basin
    # Re-organize reaches by parent
    # Remove reaches where all confluencing are outside basin(s) of interest
    relevant_tree = {}
    for parent in parent_dict:
        if len(parent_dict[parent]) < 2:
            continue
        else:
            all_basins = set(full_tree[r]['basin'] for r in parent_dict[parent])
            relevant_basins = set(BASIN_IDS)
            if all_basins.intersection(relevant_basins):
                for reach in parent_dict[parent]:
                    relevant_tree[reach] = full_tree[reach]
                relevant_tree[parent] = full_tree[parent]
    print(f'removed {len(full_tree)-len(relevant_tree)} irrelevant reaches')
    return relevant_tree


def append_elevation_data(reach_dict):
    wrk_path = os.path.join(DATA_PATH, 'Stream_Points', REACH_TYPE)

    # Organize reaches by basin
    basin_reach_dict = {}
    for reach in reach_dict:
        if not reach_dict[reach]['basin'] in basin_reach_dict:
            basin_reach_dict[reach_dict[reach]['basin']] = [reach]
        else:
            basin_reach_dict[reach_dict[reach]['basin']].append(reach)

    # Find elevation by reach
    reach_elevation_dict = {}
    for basin in basin_reach_dict:
        point_path = os.path.join(wrk_path, 'HUC12_' + basin[-4:], 'stream_points_DEMvalues.shp')
        points = gpd.read_file(point_path)
        points = points.rename(columns={'RASTERVALU': 'ELEV'})  # No good reason for this.  Keeping convention with the rest of the codebase
        for reach in basin_reach_dict[basin]:
            subset = points[points['Code'] == float(reach)]
            reach_dict[reach]['min_el'] = subset['ELEV'].min()
            reach_dict[reach]['max_el'] = subset['ELEV'].max()


def append_reach_depths(reach_dict, percentiles, intervals):
    wrk_path = os.path.join(DATA_PATH, 'Output_Logbooks', 'static', REACH_TYPE, 'thresh_' + THRESHOLD)

    # Organize reaches by basin
    basin_reach_dict = {}
    for reach in reach_dict:
        if not reach_dict[reach]['basin'] in basin_reach_dict:
            basin_reach_dict[reach_dict[reach]['basin']] = [reach]
        else:
            basin_reach_dict[reach_dict[reach]['basin']].append(reach)

    # Find depths by basin
    reach_depth_dict = {}
    for basin in basin_reach_dict:
        point_path = os.path.join(wrk_path, 'HUC12_' + basin[-4:], 'RIstage_logbook.csv')
        points_csv = pd.read_csv(point_path)
        for reach in basin_reach_dict[basin]:
            reach_dict[reach]['depths'] = {}
            for percentile in percentiles:
                if percentile not in reach_dict[reach]['depths']:
                    reach_dict[reach]['depths'][percentile] = {}
                for ri in intervals:
                    stage = points_csv['STAGE'][(points_csv['RI'] == ri) &
                                                (points_csv['REACH'] == float(reach)) &
                                                (points_csv['RESOLUTION'] == float(RESOLUTION))]
                    reach_dict[reach]['depths'][percentile][ri] = np.percentile(stage, 100 - percentile)


def calculate_backwater_elevations(reach_dict):
    # Re-organize reaches by parent
    reaches_by_parent = {}
    for reach in reach_dict:
        if reach_dict[reach]['parent'] not in reaches_by_parent:
            reaches_by_parent[reach_dict[reach]['parent']] = [reach]
        else:
            reaches_by_parent[reach_dict[reach]['parent']].append(reach)
    # remove direct connections
    direct_connects = [p for p in reaches_by_parent if len(reaches_by_parent[p]) < 2]
    for parent in direct_connects:
        reaches_by_parent.pop(parent, None)

    # segregate confluencing reaches by where they enter the parent.  Then calculate backwater depths
    for parent in reaches_by_parent:
        # Segregate
        top_confluence = []
        mid_confluence = []
        min_el = reach_dict[parent]['min_el']
        max_el = reach_dict[parent]['max_el']
        for reach in reaches_by_parent[parent]:
            tmp_min = reach_dict[reach]['min_el']
            if abs(tmp_min - max_el) <= 0.1:
                top_confluence.append(reach)
            elif min_el <= tmp_min <= max_el:
                mid_confluence.append(reach)

        # calculate backwater elevations for top-confluencing
        if len(top_confluence) > 1:
            max_depths = {}
            for reach in top_confluence:
                for percentile in reach_dict[reach]['depths']:
                    if percentile not in max_depths:
                        max_depths[percentile] = {}
                    for ri in reach_dict[reach]['depths'][percentile]:
                        tmp_depth = reach_dict[reach]['depths'][percentile][ri]
                        if ri not in max_depths[percentile]:
                            max_depths[percentile][ri] = tmp_depth
                        elif tmp_depth > max_depths[percentile][ri]:
                            max_depths[percentile][ri] = tmp_depth
            for reach in top_confluence:
                if 'bw_depths' not in reach_dict[reach]:
                    reach_dict[reach]['bw_depths'] = {}
                    reach_dict[reach]['bw_elevations'] = {}
                for percentile in max_depths:
                    if percentile not in reach_dict[reach]['bw_depths']:
                        reach_dict[reach]['bw_depths'][percentile] = {}
                        reach_dict[reach]['bw_elevations'][percentile] = {}
                    for ri in max_depths[percentile]:
                        reach_dict[reach]['bw_depths'][percentile][ri] = max_depths[percentile][ri]
                        reach_dict[reach]['bw_elevations'][percentile][ri] = max_depths[percentile][ri] + reach_dict[reach]['min_el']

        # calculate backwater elevations for mid-confluencing
        max_depths = reach_dict[parent]['depths']
        for reach in mid_confluence:
            if 'bw_depths' not in reach_dict[reach]:
                reach_dict[reach]['bw_depths'] = {}
                reach_dict[reach]['bw_elevations'] = {}
            for percentile in max_depths:
                if percentile not in reach_dict[reach]['bw_depths']:
                    reach_dict[reach]['bw_depths'][percentile] = {}
                    reach_dict[reach]['bw_elevations'][percentile] = {}
                for ri in max_depths[percentile]:
                    reach_dict[reach]['bw_depths'][percentile][ri] = max_depths[percentile][ri]
                    reach_dict[reach]['bw_elevations'][percentile][ri] = max_depths[percentile][ri] + reach_dict[reach][
                        'min_el']


def get_reach_data():
    # Parse HAND data for appropriate inundation elevations
    reaches = get_reach_tree(TMP_HARDCODED_TREE_PATH)
    append_elevation_data(reaches)
    append_reach_depths(reaches, PERCENTILES, INTERVAL_LIST)
    calculate_backwater_elevations(reaches)
    return reaches


def reformat_confluences(reaches):
    # Group reaches by basin.  Remove extraneous data.  Remove extraneous basins.

    reformatted = {}
    for reach in reaches:
        # Add reach to appropriate HUC12
        basin = reaches[reach]['basin']
        if basin not in BASIN_IDS:
            continue
        if 'bw_elevations' not in reaches[reach]:
            continue
        if basin not in reformatted:
            reformatted[basin] = {}

        # Reformat reach data and remove extraneous data
        for percentile in reaches[reach]['depths']:
            if percentile not in reformatted[basin]:
                reformatted[basin][percentile] = {}
            for ri in reaches[reach]['depths'][percentile]:
                if ri not in reformatted[basin][percentile]:
                    reformatted[basin][percentile][ri] = {}
                reformatted[basin][percentile][ri][reach] = {
                    'elevation': reaches[reach]['bw_elevations'][percentile][ri],
                    'depth': reaches[reach]['bw_depths'][percentile][ri]}

    return reformatted


def process_by_basin(reach_dict):
    for basin in reach_dict:
        print(f'processing basin {basin}')
        process_basin(basin, reach_dict[basin])


def process_basin(basin, data):
    #dem_path = os.path.join(DATA_PATH, 'DEM_derivatives', f'thresh_{THRESHOLD}', 'HUC12_' + basin[-4:], 'dem_filled.tif')
    dem_path = os.path.join(DATA_PATH, 'DEM_derivatives', f'thresh_{THRESHOLD}', 'HUC12_' + basin[-4:],
                            'dem_filled.tif')
    thiessen_path = os.path.join(DATA_PATH, 'Thiessen_Polygons', 'NHD', 'HUC12_' + basin[-4:], f'thiessen_{str(RESOLUTION)}m.tif')
    hand_path = os.path.join(DATA_PATH, 'HAND', f'{basin}_thresh-{THRESHOLD}_HAND_imputed{str(RESOLUTION)}m.tif')

    dem = Raster(dem_path)
    thiessen = Raster(thiessen_path)
    hand = Raster(hand_path)

    for percentile in data:
        print(f'--processing {percentile} percentile')
        for ri in data[percentile]:
            print(f'---processing {ri}-yr inundation')
            inundation_path = os.path.join(DATA_PATH, 'Inundation_Rasters', 'static', REACH_TYPE, f'thresh_{THRESHOLD}',
                                           f'HUC12_{basin[-4:]}',
                                           f'Q{str(ri)}_Depth_{str(RESOLUTION)}m_p{str(percentile)}.tif')
            inundation_raster = Raster(inundation_path)
            parameters = {
                'cols': inundation_raster.cols,
                'rows': inundation_raster.rows,
                'pixel_height': inundation_raster.pixel_height,
                'pixel_width': inundation_raster.pixel_width,
                'origin_x': inundation_raster.origin_x,
                'origin_y': inundation_raster.origin_y,
                'crs': inundation_raster.crs,
                'no_data': inundation_raster.band.GetNoDataValue()
            }
            counter = 1
            for reach in data[percentile][ri]:
                print(f'----getting inundated cells for reach {reach} | {counter} / {len(data[percentile][ri])} |')
                counter += 1

                backwater_el = data[percentile][ri][reach]['elevation']
                backwater_depth = data[percentile][ri][reach]['depth']
                #mask = (thiessen.data == int(reach)) & (dem.data < backwater_el) & (hand.data < backwater_depth)
                mask = (thiessen.data == int(reach)) & (dem.data < backwater_el) & (hand.data < backwater_depth) \
                       & (dem.data != dem.band.GetNoDataValue()) & (hand.data != hand.band.GetNoDataValue())
                depth_array = backwater_el - dem.data[mask]
                base_depth_array = inundation_raster.data[mask]
                corrected_depth_array = np.maximum(depth_array, base_depth_array)

                inundation_raster.data[mask] = corrected_depth_array

            export(basin, percentile, ri, parameters, inundation_raster.data)


def export(basin, percentile, ri, parameters, data):

    wrk_path = os.path.join(DATA_PATH, 'Post-processing')
    if not os.path.exists(wrk_path):
        os.mkdir(wrk_path)
    wrk_path = os.path.join(wrk_path, 'backwaters')
    if not os.path.exists(wrk_path):
        os.mkdir(wrk_path)
    wrk_path = os.path.join(wrk_path, basin)
    if not os.path.exists(wrk_path):
        os.mkdir(wrk_path)
    out_path = os.path.join(wrk_path, f'Q{str(ri)}_Depth_{str(RESOLUTION)}m_p{str(percentile)}.tif')

    driver = gdal.GetDriverByName('GTiff')

    out_raster = driver.Create(out_path, parameters['cols'], parameters['rows'], 1, gdal.GDT_Float32)
    out_raster.SetGeoTransform((parameters['origin_x'], parameters['pixel_width'], 0, parameters['origin_y'], 0, parameters['pixel_height']))
    out_band = out_raster.GetRasterBand(1)
    out_band.SetNoDataValue(parameters['no_data'])
    out_band.WriteArray(data)
    out_raster.SetProjection(parameters['crs'])
    out_band.FlushCache()


def export_og(basin, data, parameters):
    wrk_path = os.path.join(DATA_PATH, 'Post-processing')
    if not os.path.exists(wrk_path):
        os.mkdir(wrk_path)
    wrk_path = os.path.join(wrk_path, 'backwaters')
    if not os.path.exists(wrk_path):
        os.mkdir(wrk_path)
    wrk_path = os.path.join(wrk_path, basin)
    if not os.path.exists(wrk_path):
        os.mkdir(wrk_path)

    driver = gdal.GetDriverByName('GTiff')

    for percentile in data:
        for ri in data[percentile]:
            out_path = os.path.join(wrk_path, f'Q{str(ri)}_Depth_{str(RESOLUTION)}m_p{str(percentile)}.tif')
            out_raster = driver.Create(out_path, parameters['cols'], parameters['rows'], 1, gdal.GDT_Float32)
            out_raster.SetGeoTransform((parameters['origin_x'], parameters['pixel_width'], 0, parameters['origin_y'], 0, parameters['pixel_height']))
            out_band = out_raster.GetRasterBand(1)
            out_band.WriteArray(data[percentile][ri])
            out_raster.SetProjection(parameters['crs'])
            out_band.FlushCache()


def main():
    print('Getting backwater data at confluences')
    reaches = get_reach_data()
    reaches_by_basin = reformat_confluences(reaches)

    print('post-processing rasters')
    process_by_basin(reaches_by_basin)


if __name__ == '__main__':
    main()
