"""
Title: io_helper.py
Author: Matus Tomko
Mail: matus.tomko __at__ fmph.uniba.sk
"""
import gzip
import json
import numpy as np
import multiprocessing
import os
import pickle


def prepare_recording_vector(rec_vec, i):
    """
    Prepares a recording vectors for saving.

    Parameters
    ----------
    rec_vec : neuron.hoc.HocObject
        the recording vector

    Returns
    -------
    vec : dict
        the recording vector as a dictionary
    """
    vec = {
        'section': rec_vec.section,
        'segment_x': rec_vec.segment_x,
        'vector': np.array(rec_vec.vector)
    }
    return vec


class IOHelper:
    """
    A class used to load and save data

    ...

    Attributes
    ----------
    path_saving : str
        the path to the directory where the data is stored or will be saved
    path_settings : str
        the path to the synapses .json file
    npool : str, optional
        the number of pool processes (default multiprocessing.cpu_count() - 1)

    Methods
    -------
    save_recordings(synapses, tw_vec, v_soma_vec, t_vec, dend_vecs, p_vec, d_vec, alpha_scount_vec, ta_vec,
                        apc_vec, cai_vecs, cal_ica_vecs, ina_vecs, nmda_ica_vecs)
        Saves the recorded data
    prepare_dict_recording_vectors(vecs)
        Prepares a dictionary of recording vectors for saving.
    load_synapses()
        Loads setting from the synapses.json file.
    save_synapses(synapses)
        Saves synapses to a .json file.
    load_setting()
        Loads setting from the setting.json file.
    save_setting(setting)
        Saves setting to the setting.json file.
    """

    def __init__(self, path_saving, path_settings):
        """
        Parameters
        ----------
        path_saving : str
            the path to the directory where the data will be saved
        path_settings : str
            the path to the directory with settings
        """

        self.path_saving = path_saving
        self.path_settings = path_settings
        self.npool = multiprocessing.cpu_count() - 1

        try:
            if not os.path.exists(self.path_saving):
                os.makedirs(self.path_saving)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

    def save_recordings(self, synapses, tw_vec, v_soma_vec, t_vec, dend_vecs, p_vec, d_vec, alpha_scount_vec, ta_vec,
                        apc_vec, cai_vecs, cal_ica_vecs, ina_vecs, nmda_ica_vecs):
        """
        Saves the recorded data to a dictionary structures in binary files.

        Parameters
        ----------
        synapses : dict
            the dictionary of synapses
        tw_vec : neuron.hoc.HocObject
            the time vector for synaptic weights
        v_soma_vec : neuron.hoc.HocObject
            the somatic voltage vector
        t_vec : neuron.hoc.HocObject
            the time vector for voltage
        dend_vecs : dict
            the dictionary containing voltage vectors from dendrites
        p_vec : neuron.hoc.HocObject
            the potentiation amplitude vector
        d_vec : neuron.hoc.HocObject
            the depression amplitude vector
        alpha_scount_vec : neuron.hoc.HocObject
            the integrated spike count vector
        ta_vec : neuron.hoc.HocObject
            the time vector for amplitudes
        apc_vec : neuron.hoc.HocObject
            the vector of times of fired action potentials
        cai_vecs : dict
            the dictionary containing intracellular calcium concentration vectors
        cal_ica_vecs : dict
            the dictionary containing CaL channel-mediated calcium current vectors
        ina_vecs : dict
            the dictionary containing sodium current vectors
        nmda_ica_vecs : dict
            the dictionary of NMDAR channel-mediated calcium current vectors
        """
        # a dictionary of synapses
        synapses_dict = {}
        for sec in synapses:
            synapses_list = []
            for syn in synapses[sec]:
                s = {
                    'name': str(syn.synapse),
                    'section': str(syn.section),
                    'segment_x': syn.segment_x,
                    'distance': syn.distance,
                    'weight': np.array(syn.weight_vec),
                    'input_spikes': np.array(syn.input_spikes),
                    'stimulated': syn.stimulated,
                    'receptor': syn.receptor,
                    'pathway': syn.pathway
                }
                synapses_list.append(s)
            synapses_dict[sec] = synapses_list

        print('Saving recordings...')

        # saving of synapses
        weights = {}
        weights['T'] = np.array(tw_vec)
        weights['synapses'] = synapses_dict
        pickle.dump(weights, gzip.GzipFile(self.path_saving + 'synapses.p', 'wb'))
        print('The synapses were saved in the directory: ' + self.path_saving)

        # saving of voltages
        voltages = {}
        voltages['T'] = np.array(t_vec)
        voltages['V_soma'] = np.array(v_soma_vec)
        voltages['APs'] = np.array(apc_vec)
        voltages['V_dends'] = self.prepare_dict_recording_vectors(vecs=dend_vecs)
        pickle.dump(voltages, gzip.GzipFile(self.path_saving + 'voltages.p', 'wb'))
        print('The voltages were saved in the directory: ' + self.path_saving)

        # saving of currents
        currents = {}
        currents['T'] = np.array(t_vec)
        currents['cai'] = self.prepare_dict_recording_vectors(vecs=cai_vecs)
        currents['cal_ica'] = self.prepare_dict_recording_vectors(vecs=cal_ica_vecs)
        currents['nmda_ica'] = self.prepare_dict_recording_vectors(vecs=nmda_ica_vecs)
        currents['ina'] = self.prepare_dict_recording_vectors(vecs=ina_vecs)
        pickle.dump(currents, gzip.GzipFile(self.path_saving + 'currents.p', 'wb'))
        print('The currents were saved in the directory: ' + self.path_saving)

        # saving of BCM parameters
        bcm = {}
        bcm['P_amp'] = np.array(p_vec)
        bcm['D_amp'] = np.array(d_vec)
        bcm['alpha_scount'] = np.array(alpha_scount_vec)
        bcm['T_BCM'] = np.array(ta_vec)
        pickle.dump(bcm, gzip.GzipFile(self.path_saving + 'bcm.p', 'wb'))
        print('The BCM recordings were saved in the directory: ' + self.path_saving)

    def prepare_dict_recording_vectors(self, vecs):
        """
        Prepares a dictionary of recording vectors for saving.

        Parameters
        ----------
        vecs : dict
            the dictionary of recording vectors

        Returns
        -------
        vecs_dict : dict
            the dictionary of recording vectors
        """
        pool = multiprocessing.Pool(processes=self.npool, maxtasksperchild=1)
        vecs_dict = {}
        for sec in vecs:
            saved_nmda_ica_vecs = [pool.apply(prepare_recording_vector, args=(vecs[sec][i], i)) for i in
                                   range(len(vecs[sec]))]
            vecs_dict[sec] = saved_nmda_ica_vecs
        pool.terminate()
        pool.join()
        del pool
        return vecs_dict

    def load_synapses(self):
        """
        Loads setting from the synapses.json file.

        Returns
        --------
        synapses : dict
            the dictionary of synapses

        Raises
        ------
        FileNotFoundError
        """
        try:
            with open(self.path_settings + 'synapses.json', 'r') as f:
                synapses = json.load(f)
                return synapses
        except FileNotFoundError as fnf_error:
            raise fnf_error

    def save_synapses(self, synapses):
        """
        Saves synapses to a .json file.

        Parameters
        ----------
        synapses : dict
            the dictionary of synapses
        """
        data = {}
        for sec in synapses:
            synapses_list = []
            for s in synapses[sec]:
                w = np.array(s.weight_vec)
                syn = {
                    'name': str(s.synapse),
                    'section': str(s.section),
                    'segment_x': s.segment_x,
                    'distance': s.distance,
                    'stimulated': s.stimulated,
                    'pathway': s.pathway
                }
                if len(w) > 0:
                    syn['init_weight'] = w[0]
                    syn['final_weight'] = w[-1]
                else:
                    syn['init_weight'] = s.init_weight
                    syn['final_weight'] = None
                synapses_list.append(syn)
            data[sec] = synapses_list

        with open(self.path_saving + 'synapses.json', 'w') as f:
            json.dump(data, f, indent=4)
            print('Synapses were saved in the directory: ' + self.path_saving)

    def load_setting(self):
        """
        Loads setting from the setting.json file.

        Returns
        -------
        setting : dict
            the dictionary of setting

        Raises
        ------
        FileNotFoundError
        """
        try:
            with open(self.path_settings + 'setting.json', 'r') as f:
                setting = json.load(f)
                return setting
        except FileNotFoundError as fnf_error:
            raise fnf_error

    def save_setting(self, setting):
        """Saves setting to the setting.json file.

        Parameters
        ----------
        setting : dict
            the dictionary of setting
        """
        with open(self.path_saving + 'setting.json', 'w') as f:
            json.dump(setting, f, indent=4)
            print('Setting was saved in the directory: ' + self.path_saving)
