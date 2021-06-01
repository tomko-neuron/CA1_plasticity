"""
Title: stimulation_protocols.py
Author: Matus Tomko
Mail: matus.tomko __at__ fmph.uniba.sk
"""
import numpy as np
from neuron import h
from CA1_plasticity.model.utils import Synapse


class StimulationProtocol:
    """
    A class used to set a stimulation protocol

    ...

    Attributes
    ----------
    setting : dict
        the setting dictionary
    net_cons : list
        the list of neuron.hoc.NetCons
    ppStims : list
        the list of neuron.hoc.spGen2s
    vec_stims : list
        the list of neuron.hoc.VecStims

    Methods
    -------
    create_VecStim(t_vec, synapse)
        Creates a vector of stimulus times.
    set_Dong_sequential_stimulation(synapses)
        Sets the sequential stimulation stimulation protocol using Vecstim objects for Dong et al experiments.
    set_Pavlowsky_Alarcon_HFS(synapses)
        Sets the HFS stimulation protocol using Vecstim objects for Pavlowsky & Alarcon experiments.
    set_Pavlowsky_Alarcon_LFS(synapses)
        Sets the LFS stimulation protocol using Vecstim objects for Pavlowsky & Alarcon experiments.
    set_Pavlowsky_Alarcon_PP(synapses)
        Sets the paired-pulses stimulation protocol using Vecstim objects for Pavlowsky & Alarcon experiments.
    set_ppStim(synapses)
        Sets the paired-pulses stimulation protocol using spGen2 objects.
    set_square_pulse(synapses)
        Sets the square pulse stimulation protocol using Vecstim objects.
    set_theta_burst(synapses)
        Sets the theta burst stimulation protocol using Vecstim objects.
    """

    def __init__(self, setting):
        """
        Parameters
        ----------
        setting : dict
            the dictionary containing setting
        """
        self.setting = setting
        self.net_cons = []
        self.ppStims = []
        self.vec_stims = []

    def create_VecStim(self, t_vec, synapse, weight):
        """
        Creates a vector stream of events for given synapse.

        Parameters
        ----------
        t_vec : numpy.ndarray
            the time vector
        synapse : Synapse
            synapse
        weight : float
            the stimulus weight
        """
        vec = h.Vector(t_vec)
        vec_stim = h.VecStim()
        vec_stim.play(vec)
        self.vec_stims.append(vec_stim)
        nc = h.NetCon(vec_stim, synapse.synapse, 0, 0, weight)
        nc.record(synapse.input_spikes)
        self.net_cons.append(nc)

    def set_Dong_sequential_stimulation(self, synapses):
        """
        Sets the sequential stimulation stimulation protocol using Vecstim objects for Dong et al experiments.

        Parameters
        ----------
        synapses : dict
            the dictionary containing synapses
        """
        for sec in synapses:
            for syn in synapses[sec]:
                if np.random.rand() < self.setting['protocol']['Dong_SSt']['DONG_STIMULATED_PERC']:
                    if syn.pathway == 'SCH':
                        t_start = self.setting['protocol']['Dong_SSt']['DONG_SCH_START'] + np.random.rand()
                    elif syn.pathway == 'COM':
                        t_start = self.setting['protocol']['Dong_SSt']['DONG_COM_START'] + np.random.rand()
                    else:
                        continue
                    t_stop = t_start + self.setting['protocol']['Dong_SSt']['DONG_PULSES_NUM'] * \
                             self.setting['protocol']['Dong_SSt']['DONG_INTERPULSE_INTERVAL']
                    t_vec = np.arange(t_start, t_stop, self.setting['protocol']['Dong_SSt']['DONG_INTERPULSE_INTERVAL'])
                    self.create_VecStim(t_vec=t_vec,
                                        synapse=syn,
                                        weight=self.setting['protocol']['Dong_SSt']['DONG_WEIGHT'])
                    syn.stimulated = True
                else:
                    continue

    def set_Pavlowsky_Alarcon_HFS(self, synapses):
        """
        Sets the HFS stimulation protocol using Vecstim objects for Pavlowsky & Alarcon experiments.

        Parameters
        ----------
        synapses : dict
            the dictionary containing synapses
        """
        for sec in synapses:
            for syn in synapses[sec]:
                if np.random.rand() < self.setting['protocol']['Pavlowsky_Alarcon']['HFS_STIMULATED_PERC']:
                    t_start = self.setting['protocol']['Pavlowsky_Alarcon']['HFS_START'] + np.random.rand()
                    t_vec = np.zeros(0)
                    for i in range(self.setting['protocol']['Pavlowsky_Alarcon']['HFS_TRAINS_NUM']):
                        vec = np.arange(t_start, t_start + 1000, 10)
                        t_vec = np.concatenate((t_vec, vec), axis=0)
                        t_start = t_start + 1000 + self.setting['protocol']['Pavlowsky_Alarcon']['HFS_INTERTRAIN']
                    self.create_VecStim(t_vec=t_vec,
                                        synapse=syn,
                                        weight=self.setting['protocol']['Pavlowsky_Alarcon']['HFS_WEIGHT'])
                    syn.stimulated = True
                else:
                    continue

    def set_Pavlowsky_Alarcon_LFS(self, synapses):
        """
        Sets the LFS stimulation protocol using Vecstim objects for Pavlowsky & Alarcon experiments.

        Parameters
        ----------
        synapses : dict
            the dictionary containing synapses
        """
        for sec in synapses:
            for syn in synapses[sec]:
                if np.random.rand() < self.setting['protocol']['Pavlowsky_Alarcon']['LFS_STIMULATED_PERC']:
                    t_start = self.setting['protocol']['Pavlowsky_Alarcon']['LFS_START'] + np.random.rand()
                    t_vec = np.zeros(0)
                    vec = np.arange(t_start,
                                    t_start + self.setting['protocol']['Pavlowsky_Alarcon']['LFS_STIM_LEN'],
                                    1000)
                    t_vec = np.concatenate((t_vec, vec), axis=0)
                    self.create_VecStim(t_vec=t_vec,
                                        synapse=syn,
                                        weight=self.setting['protocol']['Pavlowsky_Alarcon']['LFS_WEIGHT'])
                    syn.stimulated = True
                else:
                    continue

    def set_Pavlowsky_Alarcon_PP(self, synapses):
        """
        Sets the paired-pulses stimulation protocol using Vecstim objects for Pavlowsky & Alarcon experiments.

        Parameters
        ----------
        synapses : dict
            the dictionary containing synapses
        """
        for sec in synapses:
            for syn in synapses[sec]:
                if np.random.rand() < self.setting['protocol']['Pavlowsky_Alarcon']['PP_STIMULATED_PERC']:
                    t_start = self.setting['protocol']['Pavlowsky_Alarcon']['PP_START'] + np.random.rand()
                    t_vec = np.zeros(0)
                    for i in range(self.setting['protocol']['Pavlowsky_Alarcon']['PP_NUM']):
                        vec = [t_start, t_start + 50]
                        t_vec = np.concatenate((t_vec, vec), axis=0)
                        t_start = t_start + 1000
                    self.create_VecStim(t_vec=t_vec,
                                        synapse=syn,
                                        weight=self.setting['protocol']['Pavlowsky_Alarcon']['PP_WEIGHT'])
                    syn.stimulated = True
                else:
                    continue

    def set_ppStim(self, synapses):
        """
        Sets the paired-pulses stimulation protocol using spGen2 objects.

        Parameters
        ----------
        synapses : dict
            the dictionary containing synapses
        """
        for sec in synapses:
            for syn in synapses[sec]:
                ppStim = h.SpGen2(0.5)
                ppStim.APinburst = self.setting['protocol']['theta_burst']['PP_STIM_APINBURST']
                ppStim.t01 = self.setting['protocol']['theta_burst']['PP_STIM_T01']
                self.ppStims.append(ppStim)
                if syn.receptor == 'AMPA':
                    nc = h.NetCon(ppStim, syn.synapse, 0, 0, self.setting['protocol']['theta_burst']['PP_WEIGHT'])
                    self.net_cons.append(nc)
                    syn.weight_vec.record(nc._ref_weight[1], self.setting['simulation']['RECORDING_STEP'])
                    nc.record(syn.input_spikes)
                elif syn.receptor == 'NMDA':
                    nc = h.NetCon(ppStim, syn.synapse, 0, 0, self.setting['protocol']['theta_burst']['PP_WEIGHT'])
                    self.net_cons.append(nc)
                    nc.record(syn.input_spikes)
                syn.stimulated = True

    def set_square_pulse(self, synapses):
        """
        Sets the square pulse stimulation protocol using Vecstim objects.

        Parameters
        ----------
        synapses : dict
            the dictionary containing synapses
        """
        for sec in synapses:
            for syn in synapses[sec]:
                t_start = self.setting['protocol']['square_pulses']['SQ_START'] + np.random.rand()
                t_vec = np.zeros(0)

                for i in range(self.setting['protocol']['square_pulses']['SQ_PULSES_NUM']):
                    pulse_vec = np.arange(t_start,
                                          self.setting['protocol']['square_pulses']['SQ_INTERSPIKE_INTERVAL'] *
                                          self.setting['protocol']['square_pulses']['SQ_STIMULI_NUM'] + t_start,
                                          self.setting['protocol']['square_pulses']['SQ_INTERSPIKE_INTERVAL'])
                    t_vec = np.concatenate((t_vec, pulse_vec), axis=0)
                    t_start = t_start + self.setting['protocol']['square_pulses']['SQ_INTERSPIKE_INTERVAL'] * \
                              self.setting['protocol']['square_pulses']['SQ_STIMULI_NUM'] + \
                              self.setting['protocol']['square_pulses']['SQ_INTERPULSE_INTERVAL']

                self.create_VecStim(t_vec=t_vec,
                                    synapse=syn,
                                    weight=self.setting['protocol']['square_pulses']['SQ_WEIGHT'])
                syn.stimulated = True

    def set_theta_burst(self, synapses):
        """
        Sets the theta burst stimulation protocol using Vecstim objects.

        Parameters
        ----------
        synapses : dict
            the dictionary containing synapses
        """
        for sec in synapses:
            for syn in synapses[sec]:
                t_start = self.setting['protocol']['theta_burst']['TB_START'] + np.random.rand()
                t_vec = np.zeros(0)
                for i in range(self.setting['protocol']['theta_burst']['TB_BURSTS_NUM']):
                    t_stop = t_start + 1 + (self.setting['protocol']['theta_burst']['TB_STIMULI_NUM'] - 1) * \
                             self.setting['protocol']['theta_burst']['TB_INTERSPIKE_INTERVAL']
                    burst_vec = np.arange(t_start,
                                          t_stop,
                                          self.setting['protocol']['theta_burst']['TB_INTERSPIKE_INTERVAL'])
                    t_vec = np.concatenate((t_vec, burst_vec), axis=0)
                    t_start = t_start + (
                            self.setting['protocol']['theta_burst']['TB_STIMULI_NUM'] - 1) * \
                              self.setting['protocol']['theta_burst']['TB_INTERSPIKE_INTERVAL'] + \
                              self.setting['protocol']['theta_burst']['TB_INTERBURST_INTERVAL']
                vec = h.Vector(t_vec)
                vec_stim = h.VecStim()
                vec_stim.play(vec)
                self.vec_stims.append(vec_stim)
                if syn.receptor == 'AMPA':
                    nc = h.NetCon(vec_stim, syn.synapse, 0, 0, self.setting['protocol']['theta_burst']['TB_WEIGHT'])
                    self.net_cons.append(nc)
                    # syn.weight_vec.record(nc._ref_weight[1], self.setting['simulation']['RECORDING_STEP'])
                    nc.record(syn.input_spikes)
                elif syn.receptor == 'NMDA':
                    nc = h.NetCon(vec_stim, syn.synapse, 0, 0, self.setting['protocol']['theta_burst']['TB_WEIGHT'])
                    self.net_cons.append(nc)
                    nc.record(syn.input_spikes)
                syn.stimulated = True
