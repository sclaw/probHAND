from osgeo import gdal, ogr
import math

FLOW_DIR_PATH = r"D:\HAND\mad_river\Data\Missisquoi_River\DEM_Derivatives\HUC8_derivatives\flowdir_buff.tif"
STREAM_PATH = r"D:\HAND\mad_river\Data\Missisquoi_River\DEM_Derivatives\HUC8_derivatives\MSQ_thresh-5179976_stream_network_5m.tif"
APPEND_PTS_PATH = r"D:\HAND\mad_river\Data\Missisquoi_River\working\flowacc_add-ins.shp"
OUTPUT_PATH = r"D:\HAND\mad_river\Data\Missisquoi_River\DEM_Derivatives\HUC8_derivatives\MSQ_thresh-5179976_stream_network_5m_update.tif"


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

    def geographic_pt_to_index(self, x, y):
        col = math.floor((x - self.origin_x) / self.pixel_width)
        row = math.floor((y - self.origin_y) / self.pixel_height)
        return row, col

    def export_raster(self, out_path):
        driver = gdal.GetDriverByName('GTiff')
        out_raster = driver.Create(out_path, self.cols, self.rows, 1, gdal.GDT_Int16)
        out_raster.SetGeoTransform((self.origin_x, self.pixel_width, 0, self.origin_y, 0, self.pixel_height))
        out_band = out_raster.GetRasterBand(1)
        out_band.WriteArray(self.data)
        out_band.SetNoDataValue(self.band.GetNoDataValue())
        out_raster.SetProjection(self.crs)
        out_band.FlushCache()

def load_points(path):
    out_list = list()

    driver = ogr.GetDriverByName("ESRI Shapefile")
    data_source = driver.Open(path, 0)
    in_layer = data_source.GetLayer()

    for feature in in_layer:
        geom = feature.GetGeometryRef()
        if geom.GetGeometryName() == 'POINT':
            for pt in geom.GetPoints():
                out_list.append(pt)
        elif geom.GetGeometryName() == 'MULTIPOINT':
            for sub_line in geom:
                for pt in sub_line.GetPoints():
                    out_list.append(pt)
        else:
            raise RuntimeError('Invalid geometry type.  Needs to be point or multipoint')
    return out_list


def taudem_dinf_to_d8(radians):
    x = math.cos(radians)
    y = math.sin(radians)
    if abs(x) > 0.3827:
        x /= abs(x)
    else:
        x = 0
    if abs(y) > 0.3827:
        y /= -abs(y)
    else:
        y = 0
    return int(y), int(x)


def main():
    init_pts = load_points(APPEND_PTS_PATH)

    directions = Raster(FLOW_DIR_PATH)
    stream = Raster(STREAM_PATH)

    counter = 1
    for pt in init_pts:
        print(f'Processing point {counter} / {len(init_pts)}')
        counter += 1

        y_ind, x_ind = directions.geographic_pt_to_index(*pt)
        while stream.data[y_ind, x_ind] != 1:
            stream.data[y_ind, x_ind] = 1
            wrk_direction = taudem_dinf_to_d8(directions.data[y_ind, x_ind])
            y_ind += wrk_direction[0]
            x_ind += wrk_direction[1]

    stream.export_raster(OUTPUT_PATH)


if __name__ == '__main__':
    main()
