import os
import itertools
import numpy as np
import pandas as pd
import statsmodels.api as sm

from scipy import interpolate

import constants, paths, utils


def main():
	### Create empty DataFrame
	c2v_df = pd.DataFrame()

	### Initialize path to data directory
	stem = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
	data_dir = os.path.join(stem, 'Data')

	### Execute contribution to variance analysis
	c2v_df = execute(c2v_df, data_dir)

	### Re-order columns in dataframe
	c2v_df = c2v_df[['WS_ID', 'RI', 'RESOLUTION', 'REACH', 
					 'CONSTANT_COEFF', 'CONSTANT_STD', 'CONSTANT_C2V', 
					 'XS_AREA_COEFF', 'XS_AREA_STD', 'XS_AREA_C2V', 
					 'MANNINGS_COEFF', 'MANNINGS_STD', 'MANNINGS_C2V', 
					 'SLOPE_COEFF', 'SLOPE_STD', 'SLOPE_C2V', 
					 'Q_COEFF', 'Q_STD', 'Q_C2V',

					 'Y_STD', 'SUM', 'R_SQUARED']]

	### Export dataframe to CSV
	csv_uri = os.path.join(data_dir, 'c2v.csv')
	c2v_df.to_csv(csv_uri, index=False)

	return None


def execute(c2v_df, data_dir):

	### Iterate through HUC 8 watersheds
	for ws in ['Winooski_River']:
	# for ws in constants.ws_id_dict.keys():
		ws_id = constants.ws_id_dict[ws]
		
		### Get HUC12 list
		logbooks_dir = os.path.join(data_dir, ws, 'Output_Logbooks')
		huc12_list = os.listdir(logbooks_dir)
		huc12_list = [huc12.split('_')[1] for huc12 in huc12_list]
		huc12_list = sorted(huc12_list)

		### Iterate through HUC 12 watersheds
		for huc12 in huc12_list:
			print ('\nAnalyzing %s_%s' %(ws_id, huc12))

			### Initialize paths to logbook CSVs
			rc_logbook_uri = os.path.join(
				logbooks_dir, 'HUC12_%s/rc_logbook.csv' %huc12)
			RIstage_logbook_uri = os.path.join(
				logbooks_dir, 'HUC12_%s/RIstage_logbook.csv'%huc12)

			### Read CSV files to Pandas dataframes
			try:
				rc_df = pd.read_csv(rc_logbook_uri)
				RIstage_df = pd.read_csv(RIstage_logbook_uri)
				huc12_missing = False
			except:
				print ('\tHUC12 is missing - skip!')
				huc12_missing = True

			if huc12_missing == False:
				itertools_obj = itertools.product(
					constants.interval_list, 
					RIstage_df['REACH'].unique(), 
					RIstage_df['RESOLUTION'].unique())

				for yr, reach, cellsize in itertools_obj:

					print ('\t%dyr RI  -  Reach %d  -  %dm Resolution' 
						%(yr, reach, cellsize))

					vars_df = get_inputvars(
						rc_df, RIstage_df, yr, reach, cellsize)

					if len(vars_df) == 0:
						print ('\t\tData missing - skip')

					else:
						c2v_dict = c2v_model(vars_df)

						### Add info to c2v_dict
						c2v_dict['WS_ID'] = '%s_%s' %(ws_id, huc12)
						c2v_dict['RI'] = yr
						c2v_dict['REACH'] = reach
						c2v_dict['RESOLUTION'] = cellsize

						### Append c2v_dict to c2v_df
						c2v_df = c2v_df.append(c2v_dict, ignore_index=True)


	return c2v_df


def get_inputvars(rc_df, RIstage_df, yr, reach, cellsize):
	### Subset RIstage_df
	RIstage_df2 = RIstage_df[(RIstage_df['RI']==yr) & 
					   (RIstage_df['REACH']==reach) & 
					   (RIstage_df['RESOLUTION']==cellsize)]

	### Get dependant variables
	stage_list = RIstage_df2['STAGE']

	if np.all(stage_list==0):
		return pd.DataFrame()

	### Get independant variables associated with RI
	Q_list = RIstage_df2['Q']
	
	### Get independant variables associated with stage
	xsarea_list = []
	mannings_list = []
	slope_list = []
	for i, s in enumerate(np.array(stage_list)):
		rc_df2 = rc_df[
			(rc_df['REACH']==reach) & 
			(rc_df['RESOLUTION']==cellsize) &
			(rc_df['N_SIM']==i)]
		if len(rc_df2) == 0:
			return pd.DataFrame()

		### Linear interpolation using closest stage values
		stage_low = float((rc_df2['STAGE'][rc_df2['STAGE'] < s]).max())
		stage_high = float((rc_df2['STAGE'][rc_df2['STAGE'] > s]).min())

		xsarea_low = float(rc_df2['XS_AREA'][rc_df2['STAGE'] == stage_low])
		xsarea_high = float(rc_df2['XS_AREA'][rc_df2['STAGE'] == stage_high])
		f = interpolate.interp1d([stage_low, stage_high], [xsarea_low, xsarea_high])
		xs_area = f(s)
		xsarea_list.append(xs_area)

		mannings_low = float(rc_df2['MANNINGS'][rc_df2['STAGE'] == stage_low])
		mannings_high = float(rc_df2['MANNINGS'][rc_df2['STAGE'] == stage_high])
		f = interpolate.interp1d([stage_low, stage_high], [mannings_low, mannings_high])
		mannings = f(s)
		mannings_list.append(mannings)

		slope_low = float(rc_df2['SLOPE'][rc_df2['STAGE'] == stage_low])
		slope_high = float(rc_df2['SLOPE'][rc_df2['STAGE'] == stage_high])
		f = interpolate.interp1d([stage_low, stage_high], [slope_low, slope_high])
		slope = f(s)
		slope_list.append(slope)

	### Store variables in Pandas dataframe
	vars_df = pd.DataFrame()
	vars_df['STAGE'] = np.array(stage_list)

	vars_df['XS_AREA'] = np.array(xsarea_list)
	vars_df['MANNINGS'] = np.array(mannings_list)
	vars_df['SLOPE'] = np.array(slope_list)
	vars_df['Q'] = np.array(Q_list)

	### Transform variables
	vars_df['XS_AREA'] = vars_df['XS_AREA'] ** 0.5
	vars_df['SLOPE'] = vars_df['SLOPE'] ** 2.0
	
	return vars_df


def c2v_model(vars_df):
	### Create dictionary to store outputs
	c2v_dict = {}
		
	### Initialize independant variables
	y = vars_df['STAGE']
	X = vars_df[['XS_AREA', 'MANNINGS', 'SLOPE', 'Q']]

	### Add constant to X data
	X = sm.add_constant(X)
	X = X.rename(columns={'const': 'CONSTANT'})

	### Fit OLS model
	ols_fit = sm.OLS(y, X).fit()

	### Get OLS outputs for each parameter
	c2v_sum = 0
	ols_coeff = ols_fit.params.to_dict()
	for x_col in X.columns:
		slope = ols_coeff[x_col]
		x_std = np.std(X[x_col], ddof=1)
		y_std = np.std(y, ddof=1)
		slope_normalized = (slope * (x_std / y_std)) 
		slope_normalized_sq = slope_normalized ** 2.0
		c2v_sum += slope_normalized_sq

		c2v_dict[x_col+'_COEFF'] = slope
		c2v_dict[x_col+'_STD'] = x_std
		c2v_dict[x_col+'_C2V'] = slope_normalized_sq

	c2v_dict['Y_STD'] = y_std
	c2v_dict['R_SQUARED'] = ols_fit.rsquared
	c2v_dict['SUM'] = c2v_sum

	return c2v_dict


if __name__ == '__main__':
	main()

