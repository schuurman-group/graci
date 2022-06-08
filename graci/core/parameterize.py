"""
The Parameterize object and its associated functions.
"""
import sys as sys
import numpy as np
import mpi4py as mpi4py
import os as os

import graci.core.libs as libs
import graci.core.driver as driver
import graci.io.chkpt as chkpt
import graci.io.parse as parse

class Parameterize:
    """Class constructor for the Parameterize object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.algorithm   = None
        self.label       = 'parameterize'

        self.hamiltonian    = ''
        self.graci_ref_data = ''
        self.target_file    = ''
        self.pthresh        = 1.e-5

        # -------------------------------
        self.target_data    = {}


    def run(self):
        """re-parameterize the Hamiltonian"""
 
        # parse the target data file
        target_data = parse_target_file()

        # parse the graci reference file
        scf_data, ci_data = parse_graci_file(target_data)

        # print the parameterize header
        # show the mappings between the target data and the data in the
        # reference checkpoint file
        output.print_param_header(target_data, ci_data)

        # set up the directory structure
        for target in target_data.keys():
            try:
                os.mkdir(target)
            except FileExistsError:
                os.rmdir(target)
                os.mkdir(target)

        # run the first pass
        eval_energies(target_data, scf_data, ci_data, run_scf=True)

        # set the initial values of the parameters using the
        # named Hamiltonian
        args = (self.hamiltonian)
        p0 = libs.lib_func('get_params', args)

        lsq_out = scipy.optimize.leastsq(err_func, p0, 
                             args = (target_data, scf_data, ci_data),
                             ftol = self.pthresh, full_output=True)
        (fitp, cov_p, info, mesg, ierr) = lsq_out  
    

        output.print_param_results(fitp, info)

        return

    #
    def err_func(params, target_data, scf_data, ci_data):
        """
        Evaluate the error function. In this case, simple RMSD

        Arguments:

        targets:      target values
        values:       current values

        Returns:
        norm of the different between target and current values
        """

        # update the Hamiltonian
        libs.lib_func('update_params', params)

        iter_ener = eval_energies(target_data, scf_data, ci_data)




    #
    def eval_energies(target_data, scf_data, ci_data, run_scf=False):
        """
        evaluate all the energies in the graci data set

        """

        # open the reference data file and get the contents
        if os.path.isFile(self.graci_ref_data):
            ref_chkpt = h5py.File(self.graci_ref_data, 'r',
                                                   libver='latest')
        else:
            print('File: '+str(self.graci_ref_data)+' not found. '+
                  ' Exiting...')
            sys.exit(1)

        # run through all the reference data objects
        for name, data in target_data.items():

            # pull the reference SCF and CI objects
            scf_name = scf_data[name]
            scf_obj  = chpt.read(scf_name, file_handle=ref_chkpt)

            ci_names = ci_data[name]
            ci_objs  = [chkpt.read(ci, file_handle=ref_chkpt) 
                                                for ci in ci_names]

            # pull the reference CI states


            # change to the appropriate sub-directory
            os.chdir(name)

            # run the scf (and generate integral files if need be)
            if run_scf:
                scf_obj.load()

            # run all the CI objects 
            for ci_obj in ci_objs:
                ci_obj.run(scf_obj)

            # pull reference states from the graci ref file 
                        


        return energies



    #
    def parse_graci_file(target_vals):
        """
        Parse the GRaCI reference data and confirm what data is present
        that can be compared to the target values

        Arguments: 
        target_vals:    dictionary containing the reference data

        Returns:
        graci_vals:     a dictionary containing the corresponding GRaCI
                        objects
        """

       # open the reference data file and get the contents
        if os.path.isFile(self.graci_ref_data):
            ref_chkpt = h5py.File(self.graci_ref_data, 'r',
                                                   libver='latest')
        else:
            print('File: '+str(self.graci_ref_data)+' not found. '+
                  ' Exiting...')
            sys.exit(1)

        # get the top-level contents of the checkpoint file
        ref_contents = chkpt.contents(file_name = self.graci_ref_data)

        ci_objs  = {}
        scf_objs = {}
        for tkey in target_vals.keys():
            ci_objs[tkey]  = []
            scf_objs[tkey] = None

            # if the molecule string is in the reference 
            # object name, add it ref_objs dict
            for name in ref_contents:
                ci_type = any([ci in name for ci in params.ci_objs])
                if name in tkey:
                    if ci_type:
                        ci_objs[tkey].append(name)
                    elif 'Scf' in name:
                        scf_objs[tkey] = name

        return scf_objs, ci_objs


    def parse_target_file():
        """
        Parse the file containing the data we're going to parameterize
        wrt to

        Arguments: None
        Returns:   A dictionary containing the target data
        """

        with open(self.target_file, 'r') as t_file:
            t_file_lines = t_file.readlines()

        target_dict = {}
        for line in tfile_lines:
            key              = line[0]
            target_dict[key] = []
            # every other value should be a floating point number
            # try to convert..
            for i in range(1,len(line)):
                try:
                    val = line[i]
                    target_dict[key].append(parse.convert_value(val))
                except:
                    print('Error parsing value as str/float: '+str(val))
                    print('line = '+str(line))
                    sys.exit(1)
        
        return target_dict

