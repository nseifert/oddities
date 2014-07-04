#!/usr/bin/env python
__author__ = 'Nathan Seifert'
__version__ = "1.0"
__license__ = "MIT"

# Splatalogue Auto-Assigner
# Written by Nathan Seifert, 2014.
# Licensed under the MIT license.

import requests
from bs4 import BeautifulSoup
import numpy as np
import re


# CONSTANT BLOCK
c_c = 299792458.0

element_set = set(['H', 'D', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                   'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
                   'V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
                   'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
                   'Cs','Ba','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn'])

def peakpick( data, threshold_min, threshold_max = 0 ):

	max_check = lambda x: data[x,1] >= data[x-1,1] and data[x,1] >= data[x+1,1]
	thres_min_check = lambda x: data[x,1] > threshold_min
	thres_max_check = lambda x: data[x,1] < threshold_max

	if not threshold_max:
		return np.array([[data[i,0],data[i,1]] for i in range(1,len(data)-1) if
			max_check(i) and thres_min_check(i)])
	else:
	 	return np.array([[data[i,0],data[i,1]] for i in range(1,len(data)-1) if
	 		max_check(i) and thres_min_check(i) and thres_max_check(i)])



def get_splatalogue_entry( freq, range, num=1 ):
    wvln_max = c_c / ( (freq - range) * 1.0E6)
    wvln_min = c_c / ( (freq + range) * 1.0E6)

    r = requests.get('http://find.nrao.edu/splata-slap/slap?REQUEST=queryData&WAVELENGTH=%.9f/%.9f'
                     %(wvln_min,wvln_max))

    entries = BeautifulSoup(r.text,"xml").find("TABLEDATA").findAll("TR")

    results = {}
    i = 0 # Just an easy key iterator
    for val in entries:

        entry = val.findAll("TD")

        if entry:
            results[str(i)] = ( str(entry[4]).replace("<TD>","").replace("</TD>",""), # name
                                str(entry[4]).split(" ")[0].replace("<TD>","").replace("</TD>",""), # formula
                                float(str(entry[3]).replace("<TD>","").replace("</TD>","")), # frequency
                                str(entry[-2]).replace("<TD>","").replace("</TD>","") ) # QNs
            i += 1

    if num > i:
        num = i

    return sorted([results[key] for key in results], key = lambda k: np.abs(k[2]-freq))[:num]

def filter_results(results_list, filter_list):

    filtered_results = []
    for entry in results_list:

        # sanitize atomic labels
        entry_atoms = re.findall('[A-Z][^A-Z]*',entry[1])
        no_parens = [re.sub(r'\([^)]*\)', '', atom) for atom in entry_atoms]
        no_nums = [''.join(i for i in atom if i.isalpha()) for atom in no_parens]


        if not set(no_nums).intersection(filter_list) and len(no_nums) > 1:
            filtered_results.append(entry)

    return filtered_results


def process_peakpick( ppick, num_results, filter_list = [],  range=1.0 ):
    results = []
    residuals = []

    print 'Beginning ppick processing.... There are %d peaks to process' %(ppick.shape[0])
    for i, entry in enumerate(ppick):
        print i

        if filter_list:
            temp = filter_results(get_splatalogue_entry( entry[0], range, num_results), filter_list)

        else:
            temp = get_splatalogue_entry( entry[0], range, num = num_results)

        results.append(temp)
        residuals.append([res[2] - entry[0] for res in temp])
    print '\n\n\n ----------- \n Done processing peakpick...'
    return results, residuals


def write( ppick, results, residuals):


    OUTPUT_STR = ""

    HEADER = "{:<15} \t".format('EXPT FREQ. (MHz)')
    RESULT_HEADER = "{:<20} \t {:<20} \t {:<10} \t {:<40} \t | ".format("MOLECULE","SPLAT FREQ (MHz)", "SPLAT-EXPT (MHz)", "QNs")
    HEADER += RESULT_HEADER * max([len(entry) for entry in results])

    OUTPUT_STR += HEADER + "\r\n"

    for i in range( ppick.shape[0] ):
        # LINE = "%20.4f \t" %ppick[i,0]
        LINE = "{:<20}".format(round(ppick[i,0],4))

        if len(results[i]) == 0:
            LINE += "NO RESULTS -- U-LINE?"
        else:
            for j in range( len(results[i]) ):
                LINE += "{:<20.20} \t {:<20} \t {:<10.4} \t {:<40.40} \t | ".format(results[i][j][0],
                                                                                    results[i][j][2],
                                                                                    residuals[i][j],
                                                                                    results[i][j][3])

        OUTPUT_STR += LINE + "\r\n"

    return OUTPUT_STR


if __name__ == "__main__":

    print 'Starting Splatalogue-powered auto assigner....'

    input_spec_filename = "C:\\ROT\\Data\\CH3CN-H2S_250k_avgs_20_GHz_FT.txt"
    input_spec = np.loadtxt(input_spec_filename)

    input_inten_thres_min = 0.01

    """Inputs for peakpick(spec, inten_min, inten_max = 0):

    spec: two-column spectrum file with columns (freq,inten)
    inten_min: float specifying your minimum intensity cutoff:
    inten_max (default: 0): if 0, then no maximum cutoff. Otherwise, will cutoff peakpick above inten_max
    """
    input_spec_peakpick = peakpick(input_spec, input_inten_thres_min)


    # This creates a set of elements without C, H, N and S.
    atom_filter_list = element_set - set(['C','H','N','S'])

    """
     Inputs for process_peakpick:
         process_peakpick( peakpick_file, num_results, filter_list=[], range=1.0)

     peakpick_file: two column peakpick file outputted from peakpick() as above
     num_results: Maximum integer number of nearest Splatalogue results returned for each expt line
     filter_list: list of atoms you DON'T want in your results (default empty list)
     range: Frequency search range for Splatalogue results, e.g. expt_freq +/- range (default 1.0)


     Returns: results -- a list of lists of all results for each entry in peakpick_file, where
     peakpick[i] <---> results[i] (where results[i] is a list of up to num_results entries)
              residuals --- list of frequency OMCs of splat_freq - expt_freq
    """
    results, residuals = process_peakpick(input_spec_peakpick[:50], 5, atom_filter_list)


    """
    write(peakpick_file, results_file, residuals_file):

    Input the peakpick from before, as well as the returned results and residuals files from
    your process_peakpick call.

    Returns: a formatted string with all results tabulated. Save to a txt file if you'd like.
    """

    print write(input_spec_peakpick[:50], results, residuals)


    """
    Additional note: if you want to search for a single frequency, quickly, use the
    get_splatalogue_entry( freq, range, num=1 ) function:

    freq: frequency of line in MHz
    range: +/- in MHz to search for results
    num = return num closest results (in frequency) to search frequency

    Returns: A list of tuples, with each tuple a single result. Sorted by ascending OMC.
    """"
    freq = 7680.0108
    single_freq_result = get_splatalogue_entry(freq, range=1.0, num=10)
    for result in single_freq_result:
        print result


