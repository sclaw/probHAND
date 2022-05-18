import pandas as pd
import os
from source import constants, paths, driver, preprocessing


def export_inundation(ws, huc12, n, resolution, percentile, reach_type):
    ws_id = constants.ws_id_dict[ws]
    huc12_id = huc12.split('_')[1]
    paths_dict = paths.gen_paths(ws, ws_id, huc12_id, reach_type)
    data_dict = preprocessing.read_data2mem(paths_dict, resolution)
    RIstage_logbook_fn = 'RIstage_logbook.csv'
    RIstage_logbook_uri = os.path.join(paths_dict['logbooks_dir'], RIstage_logbook_fn)
    RIstage_df = pd.read_csv(RIstage_logbook_uri)
    for p in percentile:
        driver.map_inun(RIstage_df, data_dict, p, resolution, paths_dict)


def main():
    h12_list = ['WIN_0501', 'WIN_0502', 'WIN_0503', 'WIN_0504']
    for h12 in h12_list:
        print(f'Processing basin {h12}')
        export_inundation('Winooski_River', h12, 1000, 1, [50, 60, 70, 80, 90, 100], 'NHD')


if __name__ == '__main__':
    main()
