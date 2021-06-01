"""
Title: figure_shower.py
Author: Matus Tomko
Mail: matus.tomko __at__ fmph.uniba.sk
"""

import gzip
import multiprocessing
import os
import pickle
import matplotlib.patches as mpatches
import numpy as np
from matplotlib import pyplot as plt


def calculate_average_weights(data, secs, key, return_dict):
    """
    Calculates average weights.

    Parameters
    ----------
    data : dict
        the dictionary containing the data
    secs : list
        the list of sections names
    key : str
        the key for a return dictionary
    return_dict : dict
        the shared dictionary
    """
    stim_weights = []
    unstim_weights = []
    for sec in data['synapses']:
        if sec in secs:
            for syn in data['synapses'][sec]:
                if syn['stimulated']:
                    stim_weights.append(syn['weight'])
                else:
                    unstim_weights.append(syn['weight'])

    if len(stim_weights) > 0:
        stim_weights = np.array(stim_weights)
        avg_stim_weights = np.mean(stim_weights, axis=0)
        return_dict[key + ' stim'] = avg_stim_weights

    if len(unstim_weights) > 0:
        unstim_weights = np.array(unstim_weights)
        avg_unstim_weights = np.mean(unstim_weights, axis=0)
        return_dict[key + ' unstim'] = avg_unstim_weights

    return_dict['T'] = data['T']


def calculate_average_SCH_COM_weights(data, secs, layer, pathway, return_dict):
    """
    Calculates average weights.

    Parameters
    ----------
    data : dict
        the dictionary containing the data
    secs : list
        the list of sections names
    layer : str
        the layer (ori, rad, lm)
    pathway : str
        the pathway (SCH, COM, LM)
    return_dict : dict
        the shared dictionary
    """
    weights = []
    for sec in secs:
        for syn in data['synapses'][sec]:
            if syn['pathway'] == pathway:
                weights.append(syn['weight'])

    weights = np.array(weights)
    av_weights = np.mean(weights, axis=0)
    if pathway is None:
        return_dict[layer] = av_weights
    else:
        return_dict[layer + '_' + pathway] = av_weights
    return_dict['T'] = data['T']


class FiguresShower:
    """
    A class used to show figures

    ...

    Attributes
    ----------
    path_recordings : str
        path to the directory containing recordings
    path_saving : str
        the path to the directory where the figures will be saved
    setting : dict
        the dictionary containing setting
    save_figures : bool
        the flag indicates whether the figures will be saved (default False)

    COM : Patch
    LM : Patch
    SCH : Patch
    ori_stim : Patch
    ori_unstim : Patch
    rad_stim : Patch
    rad_unstim : Patch
    lm_stim : Patch
    lm_unstim : Patch

    ori_secs : list
    rad_secs : list
    lm_secs : list

    Methods
    -------
    show_alpha_scount()
        Shows integrated spike count scaled by alpha.
    show_d_p_amplitude()
        Shows depression and potentiation amplitudes.
    show_cal_ica()
        Shows CaL channel-mediated calcium current traces for each section.
    show_intracellular_calcium()
        Shows intracellular calcium concentration traces for each section.
    show_na_current()
        Shows Na current traces for each section.
    show_nmda_ica_current()
        Shows NMDA channel-mediated calcium current traces for each section.
    show_average_weights(secs, keys)
        Shows average weights for given sections.
    show_average_weights_change( secs, keys, baseline)
        Shows evolution of average weights in time.
    show_IO_characteristics()
        Shows average input and output firing frequencies.
    show_SCH_COM_average_weights()
        Shows average Schaffer and commissural weights.
    show_SCH_COM_average_weights_change(baseline)
        Shows evolution of average Schaffer and commissural weights in time.
    show_SCH_COM_weights()
        Shows synaptic weights for each section labeled as Schaffer, commissural or perforant.
    show_weights()
        Shows synaptic weights for each section.
    show_weights_distance(t=0)
        Shows the percentage change for each synapse between the initial weight and the weight at time t
        as a scatter plot. X-axis represents distance of a synapse from the soma.
    show_dendritic_voltage()
        Shows dendritic voltage for each section.
    show_input_spikes()
        Shows input spikes.
    show_somatic_voltage()
        Shows the somatic voltage.
    """

    def __init__(self, path_recordings, path_saving, setting, save_figures=False):
        """
        Parameters
        ----------
        path_recordings : str
            path to the directory containing recordings
        path_saving : str
            the path to the directory where the figures will be saved
        setting : dict
            the dictionary containing setting
        save_figures : bool
            the flag indicates whether the figures will be saved (default False)
        """
        self.path_recordings = path_recordings
        self.path_saving = path_saving
        self.setting = setting
        self.save_figures = save_figures

        self.COM = mpatches.Patch(color='blue', label='Commissural')
        self.LM = mpatches.Patch(color='green', label='Perforant')
        self.SCH = mpatches.Patch(color='red', label='Schaffer')
        self.ori_stim = mpatches.Patch(color='tomato', label='ori stim')
        self.ori_unstim = mpatches.Patch(color='salmon', label='ori unstim')
        self.rad_stim = mpatches.Patch(color='navy', label='rad stim')
        self.rad_unstim = mpatches.Patch(color='mediumblue', label='rad unstim')
        self.lm_stim = mpatches.Patch(color='forestgreen', label='lm stim')
        self.lm_unstim = mpatches.Patch(color='limegreen', label='lm unstim')

        self.ori_secs = ['oriprox1', 'oriprox2', 'oridist1_1', 'oridist1_2', 'oridist2_1', 'oridist2_2']
        self.rad_secs = ['radTprox1', 'radTprox2', 'radTmed1', 'radTmed2', 'radTdist1', 'radTdist2',
                         'rad_t1', 'rad_t3', 'rad_t2']
        self.lm_secs = ['lm_thick1', 'lm_medium1', 'lm_thin1', 'lm_thick2', 'lm_medium2', 'lm_thin2']

        try:
            if not os.path.exists(self.path_saving):
                os.makedirs(self.path_saving)
        except OSError as e:
            if e.errno != 17:
                raise
            pass

        plt.rc('font', family='sans-serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')

    def show_alpha_scount(self):
        """Shows integrated spike count scaled by alpha."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'bcm.p', 'rb'))

        fig = plt.figure(figsize=[15, 4])
        plt.plot(data['T_BCM'], data['alpha_scount'])
        plt.xlabel('Time (min)')
        plt.ylabel('<c> * alpha')
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'alpha_scount' + '.svg', format='svg')
            print('The Integrated spike count scaled by alpha was saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()

    def show_d_p_amplitude(self):
        """Shows depression and potentiation amplitudes."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'bcm.p', 'rb'))

        fig = plt.figure(figsize=[15, 4])
        plt.plot(data['T_BCM'], data['D_amp'], label='Depression amplitude')
        plt.plot(data['T_BCM'], data['P_amp'], label='Potentiation amplitude')
        plt.xlabel('Time (min)')
        plt.ylabel('Amplitude')
        plt.legend(loc='upper right')
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'amplitudes' + '.svg', format='svg')
            print('The Amplitudes were saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()

    def show_cal_ica(self):
        """Shows CaL channel-mediated calcium current traces for each section."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'currents.p', 'rb'))
        vecs = data['cal_ica']
        for sec in vecs:
            fig = plt.figure(figsize=[12, 4])
            for vec in vecs[sec]:
                plt.plot(data['T'], vec['vector'], label=sec + '(' + str(vec['segment_x']) + ')')
            plt.xlabel('Time (ms)')
            plt.ylabel('I(Ca Cal) (nA)')
            lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            fig.tight_layout()
            if self.save_figures:
                plt.savefig(self.path_saving + 'cal_ica_' + sec + '.svg', format='svg',
                            bbox_extra_artists=(lgd,), bbox_inches='tight')
                print('The Ca current was saved in the directory: ' + self.path_saving)
            plt.show()
            plt.close()

    def show_intracellular_calcium(self):
        """Shows intracellular calcium concentration traces for each section."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'currents.p', 'rb'))
        vecs = data['cai']
        for sec in vecs:
            fig = plt.figure(figsize=[12, 4])
            for vec in vecs[sec]:
                plt.plot(data['T'], vec['vector'], label=sec + '(' + str(vec['segment_x']) + ')')
            plt.xlabel('Time (ms)')
            plt.ylabel('Intracellular Ca concentration (uM)')
            lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            fig.tight_layout()
            if self.save_figures:
                plt.savefig(self.path_saving + 'cai_' + sec + '.svg', format='svg',
                            bbox_extra_artists=(lgd,), bbox_inches='tight')
                print('The free calcium was saved in the directory: ' + self.path_saving)
            plt.show()
            plt.close()

    def show_na_current(self):
        """Shows Na current traces for each section."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'currents.p', 'rb'))
        vecs = data['ina']
        for sec in vecs:
            fig = plt.figure(figsize=[12, 4])
            for vec in vecs[sec]:
                plt.plot(data['T'], vec['vector'], label=sec + '(' + str(vec['segment_x']) + ')')
            plt.xlabel('Time (ms)')
            plt.ylabel('I(Na) (nA)')
            lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            fig.tight_layout()
            if self.save_figures:
                plt.savefig(self.path_saving + 'ina_' + sec + '.svg', format='svg',
                            bbox_extra_artists=(lgd,), bbox_inches='tight')
                print('The Na current traces was saved in the directory: ' + self.path_saving)
            plt.show()
            plt.close()

    def show_nmda_ica_current(self):
        """Shows NMDA channel-mediated calcium current traces for each section."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'currents.p', 'rb'))
        vecs = data['nmda_ica']
        for sec in vecs:
            fig = plt.figure(figsize=[12, 4])
            for vec in vecs[sec]:
                plt.plot(data['T'], vec['vector'], label=sec + '(' + str(vec['segment_x']) + ')')
            plt.xlabel('Time (ms)')
            plt.ylabel('I(NMDA Ca) (nA)')
            lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            fig.tight_layout()
            if self.save_figures:
                plt.savefig(self.path_saving + 'nmda_ca' + sec + '.svg', format='svg',
                            bbox_extra_artists=(lgd,), bbox_inches='tight')
                print('The NMDA Ca current was saved in the directory: ' + self.path_saving)
            plt.show()
            plt.close()

    def show_average_weights(self, secs, keys):
        """
        Shows average weights for given sections.

        Parameters
        ----------
        secs : list
            the list of section names lists
        keys : list
            the list of keys, labels
        """
        assert len(secs) == len(keys), 'The length of secs and keys must be the same'

        data = pickle.load(gzip.GzipFile(self.path_recordings + 'synapses.p', 'rb'))
        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        procs = []
        for i in range(len(secs)):
            p = multiprocessing.Process(target=calculate_average_weights, args=(data, secs[i], keys[i], return_dict))
            p.start()
            procs.append(p)

        for p in procs:
            p.join()
        del procs

        fig = plt.figure(figsize=[15, 4])
        for key in return_dict:
            if key != 'T':
                plt.plot(return_dict['T'], return_dict[key], label=key)

        lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        plt.title('Average weights')
        plt.xlabel('Time (ms)')
        plt.ylabel('Weight (uS)')
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'average_weights' + '.svg', format='svg',
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            print('The Average weights were saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()

    def show_average_weights_change(self, secs, keys, baseline):
        """
        Shows evolution of average weights in time.

        Parameters
        ----------
        secs : list
            the list of section names lists
        keys : list
            the list of keys, labels
        baseline : int
            the baseline
        """
        assert len(secs) == len(keys), 'The length of secs and keys must be the same'

        data = pickle.load(gzip.GzipFile(self.path_recordings + 'synapses.p', 'rb'))
        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        procs = []
        for i in range(len(secs)):
            p = multiprocessing.Process(target=calculate_average_weights, args=(data, secs[i], keys[i], return_dict))
            p.start()
            procs.append(p)

        for p in procs:
            p.join()
        del procs

        base_idx, = np.where(np.isclose(return_dict['T'], baseline))
        fig = plt.figure(figsize=[10, 3])
        for key in return_dict:
            if key != 'T':
                baseline_w = np.mean(return_dict[key][:base_idx[0]])
                return_dict[key] = [((w - baseline_w) / baseline_w) * 100 for w in return_dict[key]]
                if key == 'ori stim':
                    plt.plot(return_dict['T'], return_dict[key], label=key, color='tomato')
                elif key == 'ori unstim':
                    plt.plot(return_dict['T'], return_dict[key], label=key, color='salmon')
                if key == 'rad stim':
                    plt.plot(return_dict['T'], return_dict[key], label=key, color='navy')
                elif key == 'rad unstim':
                    plt.plot(return_dict['T'], return_dict[key], label=key, color='mediumblue')
                if key == 'lm stim':
                    plt.plot(return_dict['T'], return_dict[key], label=key, color='forestgreen')
                elif key == 'lm unstim':
                    plt.plot(return_dict['T'], return_dict[key], label=key, color='limegreen')
        plt.axhline(y=0, linewidth=1)
        plt.xlabel('Time (ms)')
        plt.ylabel('Weight (% change)')
        plt.legend()
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'weights_changes' + '.svg', format='svg')
            print('The change in average weights over time was saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()

    def show_input_spikes(self):
        """Shows input spikes."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'synapses.p', 'rb'))
        synapses = data['synapses']
        vecs = []

        for sec in synapses:
            if len(synapses[sec]) > 0:
                for syn in synapses[sec]:
                    vecs.append(syn['input_spikes'])

        fig = plt.figure(figsize=[15, 4])
        plt.eventplot(positions=vecs, orientation='horizontal')
        plt.xlabel('Time (ms)')
        plt.ylabel('Synapse')
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'input_spikes' + '.svg', format='svg')
            print('The input spikes were saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()

    def show_SCH_COM_average_weights(self):
        """Shows average Schaffer and commissural weights."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'synapses.p', 'rb'))
        secs = [self.ori_secs, self.ori_secs, self.rad_secs, self.rad_secs, self.lm_secs,
                self.ori_secs + self.rad_secs, self.ori_secs + self.rad_secs]
        layers = ['ori', 'ori', 'rad', 'rad', 'lm', 'ori_rad', 'ori_rad']
        pathways = ['SCH', 'COM', 'SCH', 'COM', None, 'SCH', 'COM']

        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        procs = []
        for i in range(len(secs)):
            p = multiprocessing.Process(target=calculate_average_SCH_COM_weights,
                                        args=(data, secs[i], layers[i], pathways[i], return_dict))
            p.start()
            procs.append(p)

        for p in procs:
            p.join()
        del procs

        fig = plt.figure(figsize=[15, 4])
        for key in return_dict:
            if key != 'T':
                plt.plot(return_dict['T'], return_dict[key], label=key)

        lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        plt.xlabel('Time (ms)')
        plt.ylabel('Weight (uS)')
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'average_SCH_COM_weights' + '.svg', format='svg',
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            print('The average weights were saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()

    def show_SCH_COM_average_weights_change(self, baseline):
        """
        Shows evolution of average Schaffer and commissural weights in time.

        Parameters
        ----------
        baseline : int
            the baseline
        """
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'synapses.p', 'rb'))
        secs = [self.ori_secs, self.ori_secs, self.rad_secs, self.rad_secs, self.lm_secs,
                self.ori_secs + self.rad_secs, self.ori_secs + self.rad_secs]
        layers = ['ori', 'ori', 'rad', 'rad', 'lm', 'ori_rad', 'ori_rad']
        pathways = ['SCH', 'COM', 'SCH', 'COM', None, 'SCH', 'COM']

        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        procs = []
        for i in range(len(secs)):
            p = multiprocessing.Process(target=calculate_average_SCH_COM_weights,
                                        args=(data, secs[i], layers[i], pathways[i], return_dict))
            p.start()
            procs.append(p)

        for p in procs:
            p.join()
        del procs

        base_idx, = np.where(np.isclose(return_dict['T'], baseline))

        fig = plt.figure(figsize=[15, 4])
        for key in return_dict:
            if key != 'T':
                baseline_w = np.mean(return_dict[key][:base_idx[0]])
                return_dict[key] = [((w - baseline_w) / baseline_w) * 100 for w in return_dict[key]]
                plt.plot(return_dict['T'], return_dict[key], label=key)

        plt.axhline(y=0, linewidth=1)
        plt.ylim(-25, 25)
        plt.xlabel('Time (ms)')
        plt.ylabel('Weight (% change)')
        lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'weights_changes' + '.svg', format='svg',
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            print('The change in average weights over the time was saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()

    def show_SCH_COM_weights(self):
        """Shows synaptic weights for each section labeled as Schaffer, commissural or perforant."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'synapses.p', 'rb'))
        synapses = data['synapses']
        for sec in synapses:
            if len(synapses[sec]) > 0:
                fig = plt.figure(figsize=[15, 4])
                plt.title(sec)
                plt.xlabel('Time (ms)')
                plt.ylabel('Weight (uS)')
                legend = False
                for syn in synapses[sec]:
                    if syn['pathway'] == 'SCH':
                        plt.plot(data['T'], syn['weight'], color='red')
                        legend = True
                    elif syn['pathway'] == 'COM':
                        plt.plot(data['T'], syn['weight'], color='blue')
                        legend = True
                    elif syn['pathway'] == 'LM':
                        plt.plot(data['T'], syn['weight'], color='green')
                        legend = True
                    else:
                        plt.plot(data['T'], syn['weight'])
                        legend = False
                if legend:
                    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left', handles=[self.SCH, self.COM, self.LM])
                fig.tight_layout()
                if self.save_figures:
                    plt.savefig(self.path_saving + 'weights_' + sec + '.svg', format='svg')
                    print('The weights were saved in the directory: ' + self.path_saving)
                plt.show()
                plt.close()

    def show_weights(self):
        """Shows synaptic weights for each section."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'synapses.p', 'rb'))
        synapses = data['synapses']
        for sec in synapses:
            if len(synapses[sec]) > 0:
                fig = plt.figure(figsize=[15, 4])
                plt.title(sec)
                plt.xlabel('Time (ms)')
                plt.ylabel('Weight (uS)')
                for syn in synapses[sec]:
                    if len(data['T']) == len(syn['weight']):
                        plt.plot(data['T'], syn['weight'])
                fig.tight_layout()
                if self.save_figures:
                    plt.savefig(self.path_saving + 'weights_' + sec + '.svg', format='svg')
                    print('The weights were saved in the directory: ' + self.path_saving)
                plt.show()
                plt.close()

    def show_weights_distance(self, t=0):
        """
        Shows the percentage change for each synapse between the initial weight and the weight at time t
        as a scatter plot. X-axis represents distance of a synapse from the soma.

        Parameters
        ----------
        t : int
            the time
        """
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'synapses.p', 'rb'))
        index, = np.where(np.isclose(data['T'], t))
        if len(index) == 0:
            t = 0.0
            index = 0
        else:
            index = index[0]
            t = data['T'][index]

        fig = plt.figure(figsize=[8, 5])
        for sec in data['synapses']:
            if sec in self.ori_secs:
                for syn in data['synapses'][sec]:
                    weight_change = ((syn['weight'][index] - syn['weight'][0]) / syn['weight'][0]) * 100
                    if syn['stimulated']:
                        plt.scatter(x=syn['distance'], y=weight_change, color='tomato')
                    else:
                        plt.scatter(x=syn['distance'], y=weight_change, color='salmon')
            elif sec in self.rad_secs:
                for syn in data['synapses'][sec]:
                    weight_change = ((syn['weight'][index] - syn['weight'][0]) / syn['weight'][0]) * 100
                    if syn['stimulated']:
                        plt.scatter(x=syn['distance'], y=weight_change, color='navy')
                    else:
                        plt.scatter(x=syn['distance'], y=weight_change, color='mediumblue')
            elif sec in self.lm_secs:
                for syn in data['synapses'][sec]:
                    weight_change = ((syn['weight'][index] - syn['weight'][0]) / syn['weight'][0]) * 100
                    if syn['stimulated']:
                        plt.scatter(x=syn['distance'], y=weight_change, color='forestgreen')
                    else:
                        plt.scatter(x=syn['distance'], y=weight_change, color='limegreen')

        plt.axhline(y=0, linewidth=1)
        plt.xlabel('Distance from the soma (um)')
        plt.ylabel('Weight (% change)')
        plt.title('The spatial distribution of synaptic weight changes (t = ' + str(t) + ' ms)')
        lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left',
                         handles=[self.ori_stim, self.ori_unstim, self.rad_stim, self.rad_unstim,
                                  self.lm_stim, self.lm_unstim])
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'wdist' + '.svg', format='svg',
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            print('The spatial ditribution of synaptic weight changes was saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()

    def show_dendritic_voltage(self):
        """Shows dendritic voltage for each section."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'voltages.p', 'rb'))
        d_vecs = data['V_dends']
        for sec in d_vecs:
            fig = plt.figure(figsize=[12, 4])
            plt.plot(data['T'], data['V_soma'], label='soma(0.5)')
            for d_vec in d_vecs[sec]:
                plt.plot(data['T'], d_vec['vector'], label=sec + '(' + str(d_vec['segment_x']) + ')')
            plt.ylim((-80, 40))
            plt.xlabel('Time (ms)')
            plt.ylabel('Voltage (mV)')
            lgd = plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
            fig.tight_layout()
            if self.save_figures:
                plt.savefig(self.path_saving + 'dendritic_voltage_' + sec + '.svg', format='svg',
                            bbox_extra_artists=(lgd,), bbox_inches='tight')
                print('The dendritic voltage was saved in the directory: ' + self.path_saving)
            plt.show()
            plt.close()

    def show_IO_characteristics(self):
        """Shows average input and output firing frequencies."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'voltages.p', 'rb'))

        apc_n = data['APs']
        freq = len(apc_n) / (self.setting['simulation']['TSTOP'] / 1000.0)

        print("Average input frequency: " + str(self.setting['netstim']['NETSTIM_FREQUENCY']))
        print("Number of output spikes: " + str(len(apc_n)))
        print("Output firing frequency: " + str(freq))
        print('-------------------------------------------------')

    def show_somatic_voltage(self):
        """Shows the somatic voltage."""
        data = pickle.load(gzip.GzipFile(self.path_recordings + 'voltages.p', 'rb'))

        fig = plt.figure(figsize=[10, 2.5])
        plt.plot(data['T'], data['V_soma'], label='soma(0.5)')
        plt.ylim((-80, 40))
        plt.xlabel('Time (min)')
        plt.ylabel('Voltage (mV)')
        plt.legend()
        fig.tight_layout()
        if self.save_figures:
            plt.savefig(self.path_saving + 'somatic_voltage' + '.svg', format='svg')
            print('The somatic voltage was saved in the directory: ' + self.path_saving)
        plt.show()
        plt.close()
