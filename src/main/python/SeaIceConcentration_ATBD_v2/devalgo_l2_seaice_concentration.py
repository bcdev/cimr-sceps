#!/usr/bin/env python3
"""
   Python script to run CIMR L2 SIC algorithm notebook from the command-line, using papermill
"""
import os
import sys
import papermill as pm
from datetime import date
from dateutil.relativedelta import *

default_l2_dir = '.'
default_l2_grid = 'ease2-2.5km-arctic'

def run_notebook(l1b_path, l2_dir = default_l2_dir, l2_grid = default_l2_grid):
    
    notebook = 'CIMR_L2_Sea_Ice_Concentration_algorithm.ipynb'
    notebook_out = notebook.replace('.ipynb','_out.ipynb')

    notebook_par = {'l1b_path': l1b_path, 'l2_dir': l2_dir, 'l2_grid': l2_grid,}

    print("Call {} with params\n\t{}".format(notebook, notebook_par))

    _ = pm.execute_notebook(notebook,notebook_out,parameters=notebook_par)
    # _ = pm.execute_notebook(notebook,notebook_out,progress_bar=True, log_output=True, parameters=notebook_par)


if __name__ == '__main__':
    import argparse

    # prepare and parse parameters to the script
    parser = argparse.ArgumentParser(prog='devalgo_l2_seaice_concentration.py',
                                     description='run CIMR L2 SIC algorithm prepared in CIMR DEVALGO')
    parser.add_argument('-i', help='Path to the CIMR L1B input file.')
    parser.add_argument('-o', help='Directory where to write the L2 SIC file', default=default_l2_dir)
    parser.add_argument('-g', help='Target grid of the L2 SIC algorithm.', default=default_l2_grid)
    args = parser.parse_args()

    #sys.path.insert(0, os.path.abspath('./') + '/Tools/')

    # run the L2 SIC algorithm via the notebook
    try:
        run_notebook(args.i, args.o, args.g)
    except pm.exceptions.PapermillExecutionError as pme:
        print(pme)
        sys.exit("Failed with Papermill Execution Error")
    except Exception as ex:
        print(ex)
        sys.exit("Failed with general exception.")

