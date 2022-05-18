import time
from scipy.spatial import Voronoi
from osgeo import ogr, gdal
import os
import numpy as np

# Hardcoded values if script is run standalone
#SRC_PATH = r"D:\HAND\mad_river\Data\Winooski_River\Stream_Polylines\NHD"
SRC_PATH = r"D:\HAND\mad_river\Data\Winooski_River\Stream_Points\NHD"
DEST_PATH = r"D:\HAND\mad_river\Data\Winooski_River\Thiessen_Polygons\NHD"
CLIP_PATH = r"D:\HAND\mad_river\Data\Winooski_River\Watershed"


class Thiessen:

    def __init__(self):
        self.unique_field = None
        self.unique_field_type = None
        self.pixel_size = None
        self.crs = None

    def extract_stream_points_from_polylines(self, path):
        out_points = []

        # Load input shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        data_source = driver.Open(path, 0)
        layer = data_source.GetLayer()

        # Get unique key metadata
        if not self.unique_field_type:
            layer_defn = layer.GetLayerDefn()
            for i in range(layer_defn.GetFieldCount()):
                field_name = layer_defn.GetFieldDefn(i).GetName()
                if field_name == self.unique_field:
                    self.unique_field_type = layer_defn.GetFieldDefn(i).GetType()
        # Get input CRS
        self.crs = layer.GetSpatialRef()

        # Extract points from lines
        for feature in layer:
            geom = feature.GetGeometryRef()
            id_code = feature.GetField(self.unique_field)
            if geom.GetGeometryName() == 'LINESTRING':
                for pt in geom.GetPoints():
                    out_points.append({'code': id_code, 'point': pt})
            elif geom.GetGeometryName() == 'MULTILINESTRING':
                for sub_line in geom:
                    for pt in sub_line.GetPoints():
                        out_points.append({'code': id_code, 'point': pt})

        return out_points

    def extract_stream_points_from_points(self, path):
        out_points = []

        # Load input shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        data_source = driver.Open(path, 0)
        layer = data_source.GetLayer()

        # Get unique key metadata
        if not self.unique_field_type:
            layer_defn = layer.GetLayerDefn()
            for i in range(layer_defn.GetFieldCount()):
                field_name = layer_defn.GetFieldDefn(i).GetName()
                if field_name == self.unique_field:
                    self.unique_field_type = layer_defn.GetFieldDefn(i).GetType()
        # Get input CRS
        self.crs = layer.GetSpatialRef()

        # Extract points from lines
        for feature in layer:
            geom = feature.GetGeometryRef()
            id_code = feature.GetField(self.unique_field)
            if geom.GetGeometryName() == 'POINT':
                for pt in geom.GetPoints():
                    out_points.append({'code': id_code, 'point': pt})
            elif geom.GetGeometryName() == 'MULTIPOINT':
                for sub_pts in geom:
                    for pt in sub_pts.GetPoints():
                        out_points.append({'code': id_code, 'point': pt})

        return out_points

    def extract_clip_points(self, path):
        out_points = []

        # Load input shapefile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        data_source = driver.Open(path, 0)
        layer = data_source.GetLayer()

        # Extract points from lines
        for feature in layer:
            geom = feature.GetGeometryRef()
            if geom.GetGeometryName() == 'POLYGON':
                ring_count = geom.GetGeometryCount()
                if ring_count > 1:
                    raise RuntimeError('HUC-12 polygons should only have one ring')
                ring = geom.GetGeometryRef(0)
                for pt in ring.GetPoints():
                    out_points.append(pt)
            elif geom.GetGeometryName() == 'MULTIPOLYGON':
                raise RuntimeError('Multipolygons not allowed for HUC-12 layers')

        return out_points

    def process_thiessen(self, network_pts, clip_pts):
        polygons_out = []

        # Calculate Voronoi regions for all stream points
        voronoi = Voronoi([pt['point'] for pt in network_pts])

        # Calculate approximate radius for finite calculations
        x_list = [pt[0] for pt in clip_pts]
        y_list = [pt[1] for pt in clip_pts]
        radius = (((max(x_list) - min(x_list)) ** 2) + ((max(y_list) - min(y_list)) ** 2)) ** 0.5
        center = voronoi.points.mean(axis=0)

        # Get Thiessen poly for each point
        for pt in range(len(network_pts)):
            poly_vertices = [voronoi.vertices[i] for i in voronoi.regions[voronoi.point_region[pt]] if i != -1]
            infinite_indices = [i for i in voronoi.regions[voronoi.point_region[pt]] if i == -1]
            if infinite_indices:
                for i in range(len(voronoi.ridge_points)):
                    if pt in voronoi.ridge_points[i] and -1 in voronoi.ridge_vertices[i]:
                        stable_point = [i for i in voronoi.ridge_vertices[i] if i != -1]
                        pt1 = voronoi.points[voronoi.ridge_points[i][0]]
                        pt2 = voronoi.points[voronoi.ridge_points[i][1]]
                        t = pt2 - pt1  # get tangent vector
                        t /= np.linalg.norm(t)  # normalize tangent vector
                        n = np.array([-t[1], t[0]])  # normal vector
                        midpoint = np.array([pt1, pt2]).mean(axis=0)
                        direction = np.sign(np.dot(midpoint - center, n)) * n  # Check which direction to project out
                        finite_approximation = midpoint + (radius * direction)

                        poly_vertices.append(finite_approximation)
                to_sort = np.asarray(poly_vertices)
                sort_center = to_sort.mean(axis=0)
                angles = np.arctan2(to_sort[:, 1] - sort_center[1], to_sort[:, 0] - sort_center[0])
                poly_vertices = np.array(to_sort)[np.argsort(angles)]
            polygons_out.append({'code': network_pts[pt]['code'], 'poly': poly_vertices})

        return polygons_out

    def export_voronoi(self, polygons, out_directory, clip_path):
        # Create working output layer
        driver = ogr.GetDriverByName("ESRI Shapefile")
        tmp_path = os.path.join(out_directory, 'thiessen_tmp.shp')
        if os.path.exists(tmp_path):
            driver.DeleteDataSource(tmp_path)
        data_source = driver.CreateDataSource(tmp_path)
        layer = data_source.CreateLayer("voronoi", self.crs, ogr.wkbMultiPolygon)
        code_field = ogr.FieldDefn(self.unique_field, self.unique_field_type)
        layer.CreateField(code_field)

        for poly in polygons:
            feature = ogr.Feature(layer.GetLayerDefn())
            code = poly['code']
            feature.SetField(self.unique_field, code)
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for pt in poly['poly']:
                ring.AddPoint(pt[0], pt[1])
            ring.CloseRings()
            geom = ogr.Geometry(ogr.wkbPolygon)
            geom.AddGeometry(ring)
            feature.SetGeometry(geom)
            layer.CreateFeature(feature)
        feature = None

        # Create output shpapefile
        clipped_path = os.path.join(out_directory, 'thiessen.shp')
        if os.path.exists(clipped_path):
            driver.DeleteDataSource(clipped_path)
        data_source_clipped = driver.CreateDataSource(clipped_path)
        layer_clipped = data_source_clipped.CreateLayer("voronoi", self.crs, ogr.wkbMultiPolygon)

        # Load HUC-12 clipper
        clipper = driver.Open(clip_path, 0)
        clipper_layer = clipper.GetLayer()

        # Clip shapefile
        ogr.Layer.Clip(layer, clipper_layer, layer_clipped)
        layer = None
        data_source = None
        driver.DeleteDataSource(tmp_path)

        # Create output Thiessen raster
        raster_out_path = os.path.join(out_directory, f'thiessen_{self.pixel_size}m.tif')
        x_min, x_max, y_min, y_max = layer_clipped.GetExtent()
        pixel_size = self.pixel_size
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)
        target_ds = gdal.GetDriverByName('GTiff').Create(raster_out_path, x_res, y_res, 1, gdal.GDT_Int32)
        target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
        target_ds.SetProjection(self.crs.ExportToWkt())
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(-9999)

        # Rasterize clipped polygon shapefile
        gdal.RasterizeLayer(target_ds, [1], layer_clipped, options=['ATTRIBUTE=Code'])

    def generate_voronoi(self, src_path, clip_path, out_directory, unique_field, output_pixel_size, in_type='points'):
        # Initiate variables
        self.unique_field = unique_field
        self.pixel_size = output_pixel_size

        # Create Voronoi polygons'
        if in_type == 'points':
            network_pts = self.extract_stream_points_from_points(src_path)
        elif in_type == 'lines':
            network_pts = self.extract_stream_points_from_polylines(src_path)
        clip_pts = self.extract_clip_points(clip_path)
        thiessen_polygons = self.process_thiessen(network_pts, clip_pts)
        self.export_voronoi(thiessen_polygons, out_directory, clip_path)


def thiessen_polygons(points_uri, polygon_uri, out_directory, in_type='points'):
    # Static variables
    UNIQUE_FIELD = 'Code'
    OUTPUT_PIXEL_SIZE = 1

    processor = Thiessen()
    processor.generate_voronoi(points_uri, polygon_uri, out_directory, UNIQUE_FIELD, OUTPUT_PIXEL_SIZE, in_type=in_type)

    return None


def main():
    #reach_list = [p[:-8] for p in os.listdir(SRC_PATH) if p[-3:] == 'shp']
    reach_list = ['WIN_0702']
    counter = 1
    times = []
    for reach in reach_list:
        t1 = time.perf_counter()
        print('{} / {}'.format(counter, len(reach_list)))
        counter += 1
        #poly_path = os.path.join(SRC_PATH, reach + '_NHD.shp')
        poly_path = os.path.join(SRC_PATH, f'HUC12_{reach[-4:]}', 'stream_points.shp')
        clip_path = os.path.join(CLIP_PATH, reach + '.shp')
        dest_path = os.path.join(DEST_PATH, f'HUC12_{reach[-4:]}')
        #thiessen_polygons(poly_path, clip_path, dest_path, in_type='lines')
        thiessen_polygons(poly_path, clip_path, dest_path)
        times.append(time.perf_counter() - t1)
    print('average time per basin was {} seconds'.format(round(sum(times) / len(times), 2)))


if __name__ == '__main__':
    main()
