# coding: utf-8

### Filename: main.py
### Language: Python v2.7
### Created by: Jesse Gourevitch

from argparse import ArgumentParser

import constants, driver, utils

import sys

#sample call: 
    # python main.py -huc8 Winooski_River -huc12_list WIN_0502 -n 2 -reach_type SGA

def main():
    """Highest-level function. Called by user.
    
    sample calls:
        python main.py -huc8 Winooski_River -huc12_list [WIN_0501,WIN_0502,WIN_0503,WIN_0504] -n 1000 -reach_type SGA
        python main.py -huc8 Winooski_River -huc12_list [WIN_0501] -n 1000 -reach_type SGA
        python main.py -huc8 Missisquoi_River -huc12_list [MSQ_0503] -n 1000 -reach_type NHD
        python main.py -huc8 Missisquoi_River -huc12_list [MSQ_0501,MSQ_0503] -n 1000 -reach_type NHD
        python main.py -huc8 Missisquoi_River -huc12_list [MSQ_0501,MSQ_0503] -n 1000 -reach_type SGA
        python main.py -huc8 Missisquoi_River -huc12_list [MSQ_0503] -n 1000 -reach_type SGA
        python main.py -huc8 Winooski_River -huc12_list WIN_0501 -n 1000 -reach_type SGA

    Parameters:
        None

    Returns:
        None
    """

    ### Initialize argument parser
    parser = ArgumentParser()

    ### Add arguments to parser
    parser.add_argument('-huc8', dest='huc8', default=None)
    parser.add_argument('-huc12_list', dest='huc12_list', default=None)
    parser.add_argument('-n', dest='n', default=None)
    parser.add_argument('-resolutions', dest='resolutions', default=None)
    parser.add_argument('-percentiles', dest='percentiles', default=None)
    parser.add_argument('-reach_type', dest='reach_type', default='NHD')
    args = parser.parse_args()


    ### Initialize list of HUC8 watersheds to run
    huc8 = args.huc8
    if huc8 not in ['Lake_Champlain', 'Lamoille_River', 'Mettawee_River', 
                    'Missisquoi_River', 'Otter_Creek', 'Winooski_River']:
        print('ERROR: Please specify a valid HUC8 watershed!')
        sys.exit()

    ### Initialize list of HUC12 watersheds to run
    ws_id = constants.ws_id_dict[huc8]

    if args.huc12_list == None:    
        huc12_list = utils.get_huc12list(huc8, ws_id)
    
    else:
        if 'START' in args.huc12_list:
            greaterthan_huc12 = args.huc12_list[6:]
            huc12_list = utils.get_huc12list(huc8, ws_id)
            greaterthan_index = huc12_list.index(greaterthan_huc12)
            huc12_list = huc12_list[greaterthan_index:]

        else:
            huc12_list = list(map(str, args.huc12_list.strip('[]').split(',')))

    ### Initialize number of iterations in Monte Carlo simulation
    if args.n == None:
        n = 1000
    else:
        n = int(args.n)

    ### Initialize list of spatial resolutions to run
    if args.resolutions == None:    
        resolutions = [7.4, 1.0] 
        resolutions = [1.0]
    else:
        resolutions = args.resolutions

    ### Initialize list of percentiles to run
    if args.percentiles == None:    
        percentiles = [5, 10, 25, 50, 75, 90, 95]
    else:
        percentiles = args.percentiles

    ### Execute driver.py script
    driver.execute(huc8, huc12_list, n, resolutions, percentiles,args.reach_type) 

    return None


if __name__ == '__main__':
    main()