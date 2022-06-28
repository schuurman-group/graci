"""
The Parameterize object and its associated functions.
"""
import sys as sys
import numpy as np
import h5py as h5py
import mpi4py.MPI as mpi
import os as os
import graci.core.libs as libs
import graci.io.chkpt as chkpt
import graci.overlaptools.overlap as overlap

#
def eval_energies():
    """
    main
    """

    #
    comm = mpi.Comm.Get_parent()
    size = comm.Get_size()
    rank = comm.Get_rank()

    libs.lib_load('bitci')
    libs.lib_load('overlap')

    hparams = np.zeros(4, dtype=float) 
    ref_file = ''
    reference = None
    molecules = None
    ref_states = None
    scf_names = None
    ci_names = None

    hparams = comm.bcast(hparams, root=0)
    ref_file = comm.bcast(ref_file, root=0)
    reference = comm.bcast(reference, root=0)
    molecules = comm.bcast(molecules, root=0)
    ref_states = comm.bcast(ref_states, root=0)
    scf_names = comm.bcast(scf_names, root=0)
    ci_names = comm.bcast(ci_names, root=0)

    energies = {}

    # open the reference data file and get the contents
    ref_name = ref_file+'.'+str(rank)
    if os.path.isfile(ref_file):
        ref_chkpt = h5py.File(ref_name, 'r', libver='latest')
    else:
        print('File: '+str(ref_name)+' not found. Exiting...')
        sys.exit(1)

    # run through all the reference data objects
    for mindx in range(len(molecules)):

        if mindx % size != rank:
            continue

        #molecule, states in target_data.items():
        molecule = molecules[mindx]

        # pull the reference SCF and CI objects
        scf_name = scf_names[molecule]
        scf_obj  = chkpt.read(scf_name, file_handle=ref_chkpt)

        ci_name  = ci_names[molecule]
        ci_refs  = [chkpt.read(ci, file_handle=ref_chkpt)
                                    for ci in ci_name]

        # change to the appropriate sub-directory
        os.chdir(str(molecule))

        # run the scf (and generate integral files if need be)
        if reference:
            # run silent
            scf_obj.verbose = False
            scf_obj.load()

        # initialize the integrals
        if scf_obj.mol.use_df:
            type_str = 'df'
        else:
            type_str = 'exact'
        libs.lib_func('bitci_int_initialize', ['pyscf', type_str,
                    scf_obj.moint_1e, scf_obj.moint_2e_eri])

        # copy the CI objects and re-run them
        ci_runs = [ci_ref.copy() for ci_ref in ci_refs]

        # run all the CI objects 
        for ci_run in ci_runs:
            # run silent
            ci_run.verbose = False
            ci_run.update_hparam(np.asarray(hparams))
            ci_run.run(scf_obj, None)

        # deallocate the integral arrays
        libs.lib_func('bitci_int_finalize', [])

        # use overlap with ref states to identify relevant states
        ener_match = identify_states(molecule, ref_states, 
                                            ci_refs, ci_runs)

        energies[molecule] = {}
        for st, ener in ener_match.items():
            energies[molecule][st] = ener

        # move back up out of this directory
        os.chdir('../')

    comm.gather(energies, root=0)
    comm.Disconnect()

#
def identify_states(molecule, ref_states, ref_ci, new_ci):
    """
    compute overlaps to identify states
    """

    eners  = {}
    smo    = np.identity(ref_ci[0].scf.nmo, dtype=float)
    # iterate over the reference states
    for iobj in range(len(ref_ci)):
        bra_lbl = ref_ci[iobj].label
        st_dict = ref_states[molecule][bra_lbl]

        bra_st  = list(st_dict.values())
        ket_st  = list(range(new_ci[iobj].n_states()))
        Smat = overlap.overlap_st(ref_ci[iobj], new_ci[iobj],
                                bra_st, ket_st, smo, 0.95, False)

        for lbl, ref_st in st_dict.items():
            Sij = np.absolute(Smat[bra_st.index(ref_st),:])
            if np.amax(Sij) <= 0.9:
                output.print_message('MAX overlap for molecule=' +
                        str(molecule) + ' state=' + str(lbl) +
                        ' is ' + str(np.amax(Sij)))
            st = ket_st[np.argmax(Sij)]
            eners[lbl] = new_ci[iobj].energies[st]

    return eners

if __name__ == '__main__':
    eval_energies()
