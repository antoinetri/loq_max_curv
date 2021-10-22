#!/usr/bin/env python
# -*- coding: utf-8 -*-
# plotting loq_curvature
# Author: Antoine Triantafyllou 2021-05
# see paper by AC Da Silva et al., in prep., 2022

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
# import kneed see in code

def import_data(filename, sheet_name=0):
    """ function to import data from csv in df"""
    data_in = pd.read_excel(filename, sheet_name)
    return data_in


def power_law(x, a, b):
    """Function to calculate the power law with constants a and b"""
    return a*np.power(x, b)


def power_law_y(y, a, b):
    """Function to calculate y in power law"""
    return np.exp((np.log(y/a))/b)


if __name__ == '__main__':
    # specify your input dir
    filename = r'Dir_of_your_input.xlsx'

    # specify your output dir
    output_dir = r'Dir_to_your_output.csv'

    # if several sheet to iterate through, list them here
    sheets = ['cloudcal', 'easycal', 'excel', 'geoexplo']

    # specify the col name for reference conc
    x_lab, y_lab = ['_conc', '_PoD']

    # chemical elements to iterate through
    elements = ['MgO', 'Al2O3', 'SiO2', 'P2O5', 'S', 'K2O', 'CaO', 'TiO2', 'Cr', 'MnO', 'Fe2O3', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Ba']

    df_out = pd.DataFrame()
    df_out['elements'] = elements

    # iterate through sheets
    for sh in sheets:
        df_in = import_data(filename, sh)
        # create some list to stock variables
        loqmcp_out, podmcp_out, comm_out, loq25_out, med_up_25, med_up_mcp = [], [], [], [], [], []

        # iterate through elements list
        for el in elements:
            df_in2 = df_in.sort_values(str(el) + x_lab)
            # sorting and ignore nan values.
            df_in3 = df_in2.dropna(subset=[str(str(el) + x_lab), str(str(el) + y_lab)])

            xdata = df_in3[str(el) + x_lab]
            ydata = df_in3[str(el) + y_lab]

            # some conditions if empty cells are found on the way
            if ydata.empty:
                print('no value in here')
                loq_mcp = np.NaN
                pod_at_mcp = np.NaN
                loq_25 = np.NaN
                median_up = np.NaN
                median_up2 = np.NaN
                comm = "no data for this element"
                loqmcp_out.append(loq_mcp)
                podmcp_out.append(pod_at_mcp)
                comm_out.append(comm)
                loq25_out.append(loq_25)
                med_up_mcp.append(median_up)
                med_up_25.append(median_up2)
                pass

            # fitting data to power line model with least square method
            else:
                pars, cov = curve_fit(f=power_law, xdata=xdata, ydata=ydata, p0=[0, 0], bounds=(-np.inf, np.inf), maxfev=10000)
                stdevs = np.sqrt(np.diag(cov))
                res = ydata - power_law(xdata, *pars)

                a_par, b_par = pars

                # ### find the maximum curvature point coord based on power line model
                x_mod = np.linspace(start=0.0001, stop=float(df_in3[str(el) + x_lab].max()), num=200)
                y_mod = [power_law(x_temp, a_par, b_par) for x_temp in x_mod]

                # import kneed lib to help finding the maximum curvature point
                from kneed import KneeLocator
                kneedle = KneeLocator(x_mod, y_mod, S=2.0, curve="convex", direction="decreasing", interp_method='polynomial')
                print('loq_val = ' + str(round(kneedle.knee, 4)))
                # print(round(kneedle.elbow, 4))
                loq_val = round(kneedle.elbow, 4)

                # Data and results plotting
                plt.plot(xdata, ydata, "ro", label="input data")
                plt.plot(x_mod, y_mod, "b-", label="fitted power line mod")

                pod_mark = 25.0
                loq_25 = power_law_y(pod_mark, a_par, b_par)
                print('loq_25 = ' + str(loq_25))
                plt.axvline(x=loq_25, c='green', ls='--', lw=1.5, label='loq at pod25')
                plt.axvline(x=loq_val, c='black', ls='--', lw=1.5, label='loq at max curvature')
                plt.title(str(sh) + "__" + str(el))

                plt.xlabel(str(sh) + '__' + str(el) + '_conc [wt%]')
                plt.ylabel(str(sh) + '__' + str(el) +'_pod [rel. %]')
                plt.ylim([-1, 3*float(df_in3[str(el) + y_lab].max())])

                #Additional conditions: testing here if the loq value is under the min quant value, meaning the power line model is extrapoled, loq cannot be trusted
                min_val = float(df_in3[str(el) + x_lab].min())
                if loq_val < min_val:
                    print('loq is lower than the min quant value : do not trust it!')
                    loq_mcp = loq_val
                    pod_at_mcp = y_mod = power_law(loq_val, a_par, b_par)
                    comm = "below quant limit"
                    loqmcp_out.append(loq_mcp)
                    podmcp_out.append(pod_at_mcp)
                    comm_out.append(comm)
                    loq25_out.append(loq_25)

                elif loq_val == 0.0:
                    print("conc vs pod data are not distributed as a power line... no loq then")
                    loq_mcp = loq_val
                    pod_at_mcp = 0.0
                    comm = "no power line trend"
                    loqmcp_out.append(loq_mcp)
                    podmcp_out.append(pod_at_mcp)
                    comm_out.append(comm)
                    loq25_out.append(loq_25)

                elif stdevs[0] > 100.0:
                    print("covariance statistics shows a bad fitting of power line with data")
                    loq_mcp = loq_val
                    pod_at_mcp = 0.0
                    comm = "power line bad fitting"
                    loqmcp_out.append(loq_mcp)
                    podmcp_out.append(pod_at_mcp)
                    comm_out.append(comm)
                    loq25_out.append(loq_25)

                else:
                    loq_mcp = loq_val
                    pod_at_mcp = y_mod = power_law(loq_val, a_par, b_par)
                    comm = "all right"
                    loqmcp_out.append(loq_mcp)
                    podmcp_out.append(pod_at_mcp)
                    loq25_out.append(loq_25)
                    comm_out.append(comm)

                # estimate the median PoD for values that are larger than MCP
                target_col = str(el) + str(x_lab)

                rslt_df = df_in3.loc[df_in3[target_col] > loq_val]
                median_up = rslt_df[str(el) + str(y_lab)].median()
                med_up_mcp.append(median_up)

                rslt_df2 = df_in3.loc[df_in3[target_col] > loq_25]
                median_up2 = rslt_df2[str(el) + str(y_lab)].median()
                med_up_25.append(median_up2)
                print("median_up_MCP = " + str(median_up))

                # additional plotting
                plt.plot([], [], ' ', label="LOQ_MCP = " + str(round(loq_val, 1)) + " wt%")
                plt.plot([], [], ' ', label="LOQ_pod25 = " + str(round(loq_25, 1)) + " wt%")
                plt.plot([], [], ' ', label="median_pod_up_MCP = " + str(round(median_up, 1)) + " %")
                plt.plot([], [], ' ', label="median_pod_up_25 = " + str(round(median_up2, 1)) + " %")

                plt.legend(fontsize='small')
                plt.savefig(os.path.join('plots', sh, str(sh) + '_' + str(el) + '.png'))
                plt.savefig(os.path.join('plots', sh, str(sh) + '_' + str(el) + '.eps'))
                plt.show()


        # create a clean df output and export it into a csv file
        df_out[sh + '_loq_mcp'] = loqmcp_out
        df_out[sh + '_pod_mcp'] = podmcp_out
        df_out[sh + '_comm_out'] = comm_out
        df_out[sh + '_loq_at_pod25'] = loq25_out
        df_out[sh + '_med_pod_up_mcp'] = med_up_mcp
        df_out[sh + '_med_pod_up_pod25'] = med_up_25
    df_out.to_csv(output_dir)

