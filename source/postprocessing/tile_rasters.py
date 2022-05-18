from osgeo import gdal, osr
import os
import numpy as np
import math
import itertools

IN_DIR = r'D:\HAND\mad_river\Data\Missisquoi_River\Post-processing\tiled'


def raster_is_empty(path):
    data_source = gdal.Open(path)
    band = data_source.GetRasterBand(1)
    nd_value = band.GetNoDataValue()
    array = band.ReadAsArray()
    return np.all([x == nd_value for x in array])


def tile_raster(path, resolution):
    # Set up output folder
    root_dir = os.path.dirname(path)
    file_name_base = os.path.split(path)[-1][:-4]
    output_dir = os.path.join(root_dir, file_name_base)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    else:
        for f in os.listdir(output_dir):
            os.remove(os.path.join(output_dir, f))


    # Import input raster
    driver = gdal.GetDriverByName('GTiff')
    data_source = gdal.Open(path)
    band = data_source.GetRasterBand(1)
    data_type = band.DataType
    nd_value = band.GetNoDataValue()
    array = band.ReadAsArray()
    raster_cols = data_source.RasterXSize
    raster_rows = data_source.RasterYSize
    geotransform = data_source.GetGeoTransform()
    origin_x = geotransform[0]
    origin_y = geotransform[3]
    pixel_width = geotransform[1]
    pixel_height = geotransform[5]
    crs = data_source.GetProjectionRef()

    # Break array into tiles
    tile_columns = list()
    start = 0
    for col in range(math.floor(raster_cols / resolution)):
        end = start + resolution
        tile_columns.append((start, end))
        start = end
    tile_columns.append((start, raster_cols))

    tile_rows = list()
    start = 0
    for row in range(math.floor(raster_rows / resolution)):
        end = start + resolution
        tile_rows.append((start, end))
        start = end
    tile_rows.append((start, raster_rows))

    # Check tile for data and export
    exported_files = list()
    for tile in itertools.product(tile_rows, tile_columns):
        sub_array = array[tile[0][0]:tile[0][1], tile[1][0]:tile[1][1]]
        if not np.all([x == nd_value for x in sub_array]):
            out_path = os.path.join(output_dir, f'{file_name_base}_{tile[0][0]}_{tile[1][0]}.tif')
            out_raster = driver.Create(out_path, tile[1][1]-tile[1][0], tile[0][1]-tile[0][0], 1, data_type)
            tmp_x_origin = (origin_x + (tile[1][0] * pixel_width))
            tmp_y_origin = (origin_y + (tile[0][0] * pixel_height))
            out_raster.SetGeoTransform((tmp_x_origin, pixel_width, 0, tmp_y_origin, 0, pixel_height))
            out_band = out_raster.GetRasterBand(1)
            out_band.WriteArray(sub_array)
            out_crs = osr.SpatialReference()
            out_crs.ImportFromWkt(crs)
            out_raster.SetProjection(out_crs.ExportToWkt())
            out_band.FlushCache()
            exported_files.append(out_path)
            print(f'exported {len(exported_files)} tiles')

    vrt_options = gdal.BuildVRTOptions(srcNodata=nd_value, VRTNodata=nd_value)
    gdal.BuildVRT(os.path.join(output_dir, file_name_base+'.vrt'), exported_files, options=vrt_options)


def batch_tile(src_dir):
    for f in os.listdir(src_dir):
        tile_raster(os.path.join(src_dir, f), 512)


def generate_vrt(src_dir):
    tile_list = list()

    for f1 in os.listdir(src_dir):
        if os.path.isdir(os.path.join(src_dir, f1)):
            print(f'processing {f1}')
            for f2 in os.listdir(os.path.join(src_dir, f1)):
                if f2[-3:] == 'vrt':
                    tile_list.append(os.path.join(src_dir, f1, f2))

    gdal.BuildVRT(os.path.join(src_dir, 'combined.vrt'), tile_list)


if __name__ == '__main__':
    #batch_tile(IN_DIR)
    generate_vrt(IN_DIR)
