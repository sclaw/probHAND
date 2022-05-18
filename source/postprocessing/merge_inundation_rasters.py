from osgeo import gdal
import os
import numpy as np


DATA_PATH = r'D:\Academic\UVM\ffi_floodplain_mapping_fall_2021\winooski_river\Post-processing\backwaters'
OUT_PATH = r'D:\Academic\UVM\ffi_floodplain_mapping_fall_2021\winooski_river\Post-processing\merged_1-17-21'
RUN_PARAMETERS = '_Depth_1m_'


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
        self.nd_value = self.band.GetNoDataValue()

    def export(self, path):
        out_raster = self.driver.Create(path, self.cols, self.rows, 1, gdal.GDT_Int16)
        out_raster.SetGeoTransform((self.origin_x, self.pixel_width, 0, self.origin_y, 0, self.pixel_height))
        out_band = out_raster.GetRasterBand(1)
        out_band.WriteArray(self.data)
        out_band.SetNoDataValue(-9999)
        out_raster.SetProjection(self.crs)
        out_band.FlushCache()


def get_paths(huc_12, intervals, percentiles, in_dir):
    path_dict = []
    for ri in intervals:
        for percentile in percentiles:
            wrk_path = os.path.join(in_dir, huc_12, f'Q{ri}{RUN_PARAMETERS}p{percentile}.tif')
            #wrk_path = os.path.join(in_dir, f'HUC12_{huc_12[-4:]}', f'Q{ri}{RUN_PARAMETERS}p{percentile}.tif')
            assert os.path.exists(wrk_path), f'{wrk_path} not found'
            path_dict.append({'ri': ri, 'percentile': percentile, 'path': wrk_path})
    return path_dict


def process_export(path_dict, out_path, param_key):
    base = Raster(path_dict[0]['path'])
    base.data.fill(-9999)

    for layer in path_dict:
        print(f'-processing {param_key} {layer[param_key]}')
        wrk_raster = Raster(layer['path'])
        mask = (wrk_raster.data != wrk_raster.nd_value)
        base.data[mask] = layer[param_key]
        wrk_raster = None

    base.export(out_path)


def merge_rasters(type, basins, intervals, percentiles, in_dir, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    for basin in basins:
        print(f'Processing basin {basin}')
        path_dict = get_paths(basin, intervals, percentiles, in_dir)
        if type == 'interval':
            out_path = os.path.join(out_dir, f'{basin}_p{percentiles[0]}.tif')
            param_key = 'ri'
        elif type == 'percentile':
            out_path = os.path.join(out_dir, f'{basin}_Q{intervals[0]}.tif')
            param_key = 'percentile'
        process_export(path_dict, out_path, param_key)


def main():
    h12_list = ['WIN_0101', 'WIN_0102', 'WIN_0103', 'WIN_0201', 'WIN_0202', 'WIN_0203', 'WIN_0204', 'WIN_0301', 'WIN_0302', 'WIN_0401', 'WIN_0402', 'WIN_0403', 'WIN_0501', 'WIN_0502', 'WIN_0503', 'WIN_0504', 'WIN_0601', 'WIN_0602', 'WIN_0603', 'WIN_0604', 'WIN_0701', 'WIN_0702', 'WIN_0703', 'WIN_0704']
    percentiles = [50]
    intervals = [500, 200, 100, 50, 25, 10, 5, 2]
    merge_rasters('interval', h12_list, intervals, percentiles, DATA_PATH, OUT_PATH)


if __name__ == '__main__':
    main()
