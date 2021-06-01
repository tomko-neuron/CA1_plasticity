"""
Title: CA1_pyramidal_cell.py
Author: Matus Tomko
Mail: matus.tomko __at__ fmph.uniba.sk
"""
from CA1_plasticity.model.utils import RecordingVector
from CA1_plasticity.model.utils import Synapse
from neuron import h
import numpy as np


class CA1PyramidalCell:
    """
    A class represents the CA1 pyramidal cell model

    ...

    Attributes
    ----------
    setting : dict
        the setting dictionary
    CA1 : neuron.hoc.HocObject
        the CA1 pyramidal cell model
    soma : neuron.hoc.HocObject
        the HocObject representing the soma
    all : list
        the list containing all sections of the model
    apical : list
        the list containing all apical dendrites sections
    basal : list
        the list containing the basal dendrites sections
    v_vec : neuron.hoc.HocObject
        the somatic voltage vector
    t_vec : neuron.hoc.HocObject
        the time vector for somatic voltage
    t_rs_vec : neuron.hoc.HocObject
        the time vector using recording step from the setting
    cai_vecs : dict
        the dictionary containing intracellular calcium concentration vectors
    cal2_ica_vecs : dict
        the dictionary containing CaL channel-mediated calcium current vectors
    dend_vecs : dict
        the dictionary containing voltage vectors from dendrites
    ina_vecs : dict
        the dictionary containing sodium current vectors
    nmda_ica_vecs : dict
        the dictionary of NMDAR channel-mediated calcium current vectors
    apc : neuron.hoc.HocObject
        the action potential counter
    apc_vec : neuron.hoc.HocObject
        the action potential times vector
    bcm : neuron.hoc.HocObject
        the BCM mechanism
    alpha_scout_vec : neuron.hoc.HocObject
        the integrated spike count scaled by alpha vector
    d_vec : neuron.hoc.HocObject
        the depression amplitude vector
    p_vec : neuron.hoc.HocObject
        the potentiation amplitude vector
    synapses : dict
        the dictionary containing synapses
    net_cons : list
        the list containing NetCons
    net_stims : list
        the list containing NetStims
    vBoxShape : neuron.hoc.HocObject
        the box organizes a collection of graphs and command panels
    shplot : neuron.hoc.HocObject
        the Shape window

    Methods
    -------
    apply_TTX(self)
        Simulates the application of TTX as a reduction of sodium channel conductance.
    apply_nimodipine():
        Simulates the application of nimodipine as a reduction of CaL channel conductance.
    connect_spontaneous_activity(lm_syn_count, ori_syn_count, rad_syn_count, random_weights)
        Connects NetStims for spontaneous activity with synapses using NetCons.
    get_synapses(sections)
        Returns a dictionary containing synapses.
    insert_AP_counter()
        Inserts an action potential counter at the soma.
    insert_BCM()
        Inserts the BCM mechanism at the soma.
    insert_AMPA_NMDA_PP_synapses(sections)
        Inserts AMPA and NMDA synapses at each segment of the sections.
    insert_json_synapses(synapses)
        Inserts synapses from a dictionary loaded from a .json file.
    insert_random_synapses(lm_syn_count, ori_syn_count, rad_syn_count)
        Randomly inserts synapses.
    set_cai_vectors(sections)
        Sets vectors for recording of intracellular calcium concentration from sections in the list.
    set_cal2_ica_vectors(sections)
        Sets vectors for recording of CaL channel-mediated calcium current from sections in the list.
    set_dendritic_voltage_vectors(sections)
        Sets vectors for recording of voltage from sections in the list.
    set_nax_ina_vectors(sections)
        Sets vectors for recording of Na channel-mediated sodium current from sections in the list.
    set_NetStim(noise)
        Creates, sets and returns a NetStim.
    set_synapse(sec, x)
        Creates, sets and returns a synapse.
    show_synapses_PP_LTP()
        Shows a box containing the model shape plot with marked recording sites and synapses.
    """

    def __init__(self, hoc_model, path_mods, setting):
        """
        Parameters
        ----------
        hoc_model : str
            the path to model in .hoc
        path_mods : str
            the path to the directory containing the .mod files
        setting : dict
            the dictionary containing setting
        """
        self.setting = setting
        h.nrn_load_dll(path_mods + 'nrnmech.dll')
        h.xopen(hoc_model)

        self.CA1 = h.CA1_PC_Tomko()
        self.soma = self.CA1.soma[0]
        self.all = list(self.CA1.all)
        self.apical = list(self.CA1.apical)
        self.basal = list(self.CA1.basal)

        self.v_vec = h.Vector().record(self.soma(0.5)._ref_v)
        self.t_vec = h.Vector().record(h._ref_t)
        self.t_rs_vec = h.Vector().record(h._ref_t, self.setting['simulation']['RECORDING_STEP'])
        self.cai_vecs = {}
        self.cal2_ica_vecs = {}
        self.dend_vecs = {}
        self.ina_vecs = {}
        self.nmda_ica_vecs = {}

        self.apc = None
        self.apc_vec = h.Vector()

        self.bcm = None
        self.alpha_scout_vec = h.Vector()
        self.d_vec = h.Vector()
        self.p_vec = h.Vector()

        self.syn_AMPA_count = 0
        self.syn_NMDA_count = 0
        self.synapses = {}
        self.net_cons = []
        self.net_stims = []

        self.vBoxShape = None
        self.shplot = None

    def apply_TTX(self):
        """Simulates the application of TTX as a reduction of sodium channel conductance."""
        for sec in self.CA1.all:
            if h.ismembrane('nax', sec=sec):
                sec.gbar_nax = sec.gbar_nax / 2

    def apply_nimodipine(self):
        """Simulates the application of nimodipine as a reduction of CaL channel conductance."""
        for sec in self.CA1.all:
            if h.ismembrane('cal', sec=sec):
                sec.gcalbar_cal = 0

    def connect_spontaneous_activity(self, lm_syn_count=None, ori_syn_count=None, rad_syn_count=None,
                                     random_weights=False):
        """
        Connects NetStims for spontaneous activity with synapses using NetCons.

        Parameters
        ----------
        lm_syn_count : int
            the number of synapses in the str. lacunosum-moleculare dendrites (default is None)
        ori_syn_count : int
            the number of synapses in the str. oriens dendrites (default is None)
        rad_syn_count : int
            the number of synapses in the str. radiatum dendrites (default is None)
        random_weights : bool
            if True, generates random weights (default is False)
        """
        # generate random weights from the normal distribution
        if random_weights and None not in {lm_syn_count, ori_syn_count, rad_syn_count} \
                and lm_syn_count + ori_syn_count + rad_syn_count == self.syn_AMPA_count:
            loc = (self.setting['initial_weights']['ORIENS_MIN_WEIGHT']
                   + self.setting['initial_weights']['ORIENS_MAX_WEIGHT']) / 2
            std = np.std([self.setting['initial_weights']['ORIENS_MIN_WEIGHT'],
                          self.setting['initial_weights']['ORIENS_MAX_WEIGHT']])
            ori_weights = np.random.normal(loc=loc, scale=std, size=ori_syn_count)

            loc = (self.setting['initial_weights']['RADIATUM_MIN_WEIGHT']
                   + self.setting['initial_weights']['RADIATUM_MAX_WEIGHT']) / 2
            std = np.std([self.setting['initial_weights']['RADIATUM_MIN_WEIGHT'],
                          self.setting['initial_weights']['RADIATUM_MAX_WEIGHT']])
            rad_weights = np.random.normal(loc=loc, scale=std, size=rad_syn_count)

            loc = (self.setting['initial_weights']['LM_MIN_WEIGHT']
                   + self.setting['initial_weights']['LM_MAX_WEIGHT']) / 2
            std = np.std([self.setting['initial_weights']['LM_MIN_WEIGHT'],
                          self.setting['initial_weights']['LM_MAX_WEIGHT']])
            lm_weights = np.random.normal(loc=loc, scale=std, size=lm_syn_count)
            o = 0
            r = 0
            l = 0
            for sec in self.synapses:
                if sec in self.setting['section_lists']['ori_secs']:
                    for syn in self.synapses[sec]:
                        syn.init_weight = ori_weights[o]
                        o = o + 1
                elif sec in self.setting['section_lists']['rad_secs']:
                    for syn in self.synapses[sec]:
                        syn.init_weight = rad_weights[r]
                        r = r + 1
                elif sec in self.setting['section_lists']['lm_secs']:
                    for syn in self.synapses[sec]:
                        syn.init_weight = lm_weights[l]
                        l = l + 1

        # connect synapses with NetStims using NetCons
        for sec in self.synapses:
            if len(self.synapses[sec]) > 0:
                for syn in self.synapses[sec]:
                    stim = self.set_NetStim(noise=1.0)
                    nc = h.NetCon(stim, syn.synapse, 0, 0, syn.init_weight)
                    if syn.receptor == 'AMPA':
                        syn.weight_vec.record(nc._ref_weight[1], self.setting['simulation']['RECORDING_STEP'])
                    nc.record(syn.input_spikes)
                    self.net_stims.append(stim)
                    self.net_cons.append(nc)

    def get_synapses(self, sections):
        """
        Returns a dictionary containing synapses.

        Parameters
        -----------
        sections : list
            the list containing the names of sections

        Returns
        -------
        dict
            a dictionary containing synapses from sections in the list
        """
        return {sec: self.synapses[sec] for sec in sections if sec in self.synapses}

    def insert_AP_counter(self):
        """Inserts an action potential counter at the soma."""

        self.apc = h.APCount(self.soma(0.5))
        self.apc.thresh = 0
        self.apc.record(self.apc_vec)

    def insert_BCM(self):
        """Inserts the BCM mechanism at the soma."""

        self.bcm = h.BCMthreshold(self.soma(0.5))
        self.bcm.d0 = self.setting['BCM']['BCM_D0']
        self.bcm.p0 = self.setting['BCM']['BCM_P0']
        self.bcm.alpha = self.setting['BCM']['BCM_ALPHA']
        self.bcm.scount0 = self.setting['BCM']['BCM_SCOUNT0']
        self.bcm.scounttau = self.setting['BCM']['BCM_SCOUNTTAU']

        self.alpha_scout_vec.record(self.bcm._ref_alpha_scount, self.setting['simulation']['RECORDING_STEP'])
        self.d_vec.record(self.bcm._ref_d, self.setting['simulation']['RECORDING_STEP'])
        self.p_vec.record(self.bcm._ref_p, self.setting['simulation']['RECORDING_STEP'])

    def insert_AMPA_NMDA_PP_synapses(self, sections):
        """
        Inserts AMPA and NMDA synapses at each segment of the sections.

        Parameters
        ----------
        sections : list
            the list containing the sections names
        """
        h.distance(0, self.soma(0.5))
        for sec in self.CA1.all:
            if sec.hname().split('.')[1] in sections:
                synapses = []
                ica_vecs = []
                for seg in sec.allseg():
                    if seg.x not in [0.0, 1.0]:
                        dist = h.distance(seg.x, sec=sec)
                        stim = self.set_NetStim(noise=1)
                        # AMPA synapse
                        hoc_ampa_syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(sec(seg.x))
                        hoc_ampa_syn.tau1 = self.setting['synapse']['AMPA_TAU1']
                        hoc_ampa_syn.tau2 = self.setting['synapse']['AMPA_TAU2']
                        hoc_ampa_syn.e = self.setting['synapse']['AMPA_E']
                        hoc_ampa_syn.start = self.setting['synapse']['AMPA_START']
                        hoc_ampa_syn.dtau = self.setting['synapse']['AMPA_DTAU']
                        hoc_ampa_syn.ptau = self.setting['synapse']['AMPA_PTAU']
                        h.setpointer(self.bcm._ref_d, 'd', hoc_ampa_syn)
                        h.setpointer(self.bcm._ref_p, 'p', hoc_ampa_syn)
                        ampa_syn = Synapse(synapse=hoc_ampa_syn, section=sec, segment_x=seg.x, distance=dist,
                                           init_weight=self.setting['synapse']['AMPA_WEIGHT'], weight_vec=h.Vector(),
                                           input_spikes_vec=h.Vector(), receptor='AMPA')
                        synapses.append(ampa_syn)
                        self.syn_AMPA_count = self.syn_AMPA_count + 1

                        # NMDA synapse
                        hoc_nmda_syn = h.Exp2SynNMDA(sec(seg.x))
                        hoc_nmda_syn.tau1 = self.setting['synapse']['NMDA_TAU1']
                        hoc_nmda_syn.tau2 = self.setting['synapse']['NMDA_TAU2']
                        hoc_nmda_syn.e = self.setting['synapse']['NMDA_E']
                        nmda_syn = Synapse(synapse=hoc_nmda_syn, section=sec, segment_x=seg.x, distance=dist,
                                           init_weight=self.setting['synapse']['NMDA_WEIGHT'], weight_vec=h.Vector(),
                                           input_spikes_vec=h.Vector(), receptor='NMDA')
                        synapses.append(nmda_syn)

                        vec = h.Vector().record(hoc_nmda_syn._ref_ica)
                        ica_vec = RecordingVector(section=sec.hname(), segment_x=seg.x, vec=vec)
                        ica_vecs.append(ica_vec)
                        self.syn_NMDA_count = self.syn_NMDA_count + 1

                self.synapses[sec.hname().split('.')[1]] = synapses
                self.nmda_ica_vecs[sec.hname().split('.')[1]] = ica_vecs

    def insert_json_synapses(self, synapses):
        """
        Inserts synapses from a dictionary loaded from a .json file.

        Parameters
        ----------
        synapses : dict
            the dictionary containing synapses loaded from a .json file
        """
        for sec in self.CA1.all:
            syn_list = []
            for s in synapses[sec.hname().split('.')[1]]:
                hoc_syn = self.set_synapse(sec=sec, x=float(s['segment_x']))
                dist = h.distance(float(s['segment_x']), sec=sec)
                syn = Synapse(synapse=hoc_syn, section=sec, segment_x=float(s['segment_x']), distance=dist,
                              init_weight=s['init_weight'], weight_vec=h.Vector(), input_spikes_vec=h.Vector(),
                              receptor='AMPA')
                syn.pathway = s['pathway']
                syn_list.append(syn)
                self.syn_AMPA_count = self.syn_AMPA_count + 1
            self.synapses[sec.hname().split('.')[1]] = syn_list

        print('Total number of AMPA synapses: ' + str(self.syn_AMPA_count))

    def insert_random_synapses(self, lm_syn_count, ori_syn_count, rad_syn_count):
        """
        Randomly inserts synapses.

        Parameters
        ----------
        lm_syn_count : int
            the number of synapses in the str. lacunosum-moleculare dendrites
        ori_syn_count : int
            the number of synapses in the str. oriens dendrites
        rad_syn_count : int
            the number of synapses in the str. radiatum dendrites
        """
        h.distance(0, self.soma(0.5))

        # generates the list of str. oriens sections and inserts synapses
        secs = np.random.choice(self.basal, ori_syn_count, p=[0.01, 0.245, 0.245, 0.01, 0.245, 0.245])
        for sec in self.basal:
            syn_list = []
            for i in range(np.count_nonzero(secs == sec)):
                x = np.random.rand()
                hoc_syn = self.set_synapse(sec=sec, x=x)
                dist = h.distance(x, sec=sec)
                syn = Synapse(synapse=hoc_syn, section=sec, segment_x=x, distance=dist,
                              init_weight=None, weight_vec=h.Vector(), input_spikes_vec=h.Vector(), receptor='AMPA')
                syn_list.append(syn)
                self.syn_AMPA_count = self.syn_AMPA_count + 1
            self.synapses[sec.hname().split('.')[1]] = syn_list

        # generates the list of str. radiatum sections and inserts synapses
        radiatum = self.apical[2:9]
        secs = np.random.choice(radiatum, rad_syn_count, p=[0.005, 0.005, 0.06, 0.06, 0.29, 0.29, 0.29])
        for sec in radiatum:
            syn_list = []
            for i in range(np.count_nonzero(secs == sec)):
                x = np.random.rand()
                hoc_syn = self.set_synapse(sec=sec, x=x)
                h.distance(0, self.soma(0.5))
                dist = h.distance(x, sec=sec)
                syn = Synapse(synapse=hoc_syn, section=sec, segment_x=x, distance=dist,
                              init_weight=None, weight_vec=h.Vector(), input_spikes_vec=h.Vector(), receptor='AMPA')
                syn_list.append(syn)
                self.syn_AMPA_count = self.syn_AMPA_count + 1
            self.synapses[sec.hname().split('.')[1]] = syn_list

        # generates the list of str. lacunosum-moleculare sections and inserts synapses
        lm = self.apical[9:]
        secs = np.random.choice(lm, lm_syn_count)
        for sec in lm:
            syn_list = []
            for i in range(np.count_nonzero(secs == sec)):
                x = np.random.rand()
                hoc_syn = self.set_synapse(sec=sec, x=x)
                h.distance(0, self.soma(0.5))
                dist = h.distance(x, sec=sec)
                syn = Synapse(synapse=hoc_syn, section=sec, segment_x=x, distance=dist,
                              init_weight=None, weight_vec=h.Vector(), input_spikes_vec=h.Vector(), receptor='AMPA')
                syn_list.append(syn)
                self.syn_AMPA_count = self.syn_AMPA_count + 1
            self.synapses[sec.hname().split('.')[1]] = syn_list

        print('Total number of AMPA synapses: ' + str(self.syn_AMPA_count))

    def set_cai_vectors(self, sections):
        """
        Sets vectors for recording of intracellular calcium concentration from sections in the list.

        Parameters
        ----------
        sections : list
            the list containing the names of sections
        """
        for sec in self.CA1.all:
            if sec.hname().split('.')[1] in sections:
                cai_vecs = []
                for seg in sec.allseg():
                    vec = h.Vector().record(seg._ref_cai)
                    cai_vec = RecordingVector(section=sec.hname(), segment_x=seg.x, vec=vec)
                    cai_vecs.append(cai_vec)
                self.cai_vecs[sec.hname().split('.')[1]] = cai_vecs

    def set_cal2_ica_vectors(self, sections):
        """
        Sets vectors for recording of CaL channel-mediated calcium current from sections in the list.

        Parameters
        ----------
        sections : list
            the list containing the names of sections
        """
        for sec in self.CA1.all:
            if sec.hname().split('.')[1] in sections and h.ismembrane('cal', sec=sec):
                cal2_ica_vecs = []
                for seg in sec.allseg():
                    vec = h.Vector().record(seg._ref_ica_cal)
                    cal2_ica_vec = RecordingVector(section=sec.hname(), segment_x=seg.x, vec=vec)
                    cal2_ica_vecs.append(cal2_ica_vec)
                self.cal2_ica_vecs[sec.hname().split('.')[1]] = cal2_ica_vecs

    def set_dendritic_voltage_vectors(self, sections):
        """
        Sets vectors for recording of voltage from sections in the list.

        Parameters
        ----------
        sections : list
            the list containing the names of sections
        """
        for sec in self.CA1.all:
            if sec.hname().split('.')[1] in sections:
                d_vecs = []
                for seg in sec.allseg():
                    vec = h.Vector().record(seg._ref_v)
                    d_vec = RecordingVector(section=sec.hname(), segment_x=seg.x, vec=vec)
                    d_vecs.append(d_vec)
                self.dend_vecs[sec.hname().split('.')[1]] = d_vecs

    def set_nax_ina_vectors(self, sections):
        """
        Sets vectors for recording of Na channel-mediated sodium current from sections in the list.

        Parameters
        ----------
        sections : list
            the list containing the names of sections
        """
        for sec in self.CA1.all:
            if sec.hname().split('.')[1] in sections and h.ismembrane('nax', sec=sec):
                na_vecs = []
                for seg in sec.allseg():
                    vec = h.Vector().record(seg._ref_ina_nax)
                    na_vec = RecordingVector(section=sec.hname(), segment_x=seg.x, vec=vec)
                    na_vecs.append(na_vec)
                self.ina_vecs[sec.hname().split('.')[1]] = na_vecs

    def set_NetStim(self, noise=0.0):
        """
        Creates, sets and returns a NetStim.

        Parameters
        ----------
        noise : float
            the fractional randomness (default is 0.0)

        Returns
        -------
        neuron.hoc.HocObject
            NetStim
        """
        stim = h.NetStim()
        stim.seed(self.setting['netstim']['NETSTIM_SEED'])
        stim.start = np.random.rand() * self.setting['netstim']['NETSTIM_START']
        if self.setting['netstim']['NETSTIM_FREQUENCY'] == 0:
            stim.interval = 1e9
        else:
            stim.interval = 1000 / self.setting['netstim']['NETSTIM_FREQUENCY']
        stim.number = self.setting['netstim']['NETSTIM_NUMBER']
        stim.noise = noise
        return stim

    def set_synapse(self, sec, x):
        """
        Creates, sets and returns a synapse.

        Parameters
        ----------
        sec : neuron.hoc.HocObject
            the section into which synapse will be inserted
        x : neuron.hoc.HocObject
            the synapse location at the section

        Returns
        -------
        neuron.hoc.HocObject
            synapse
        """
        syn = h.Exp2SynSTDP_multNNb_globBCM_intscount_precentred(sec(x))
        syn.tau1 = self.setting['synapse']['SYN_TAU1']
        syn.tau2 = self.setting['synapse']['SYN_TAU2']
        syn.e = self.setting['synapse']['SYN_E']
        syn.start = self.setting['synapse']['SYN_START']
        syn.dtau = self.setting['synapse']['SYN_DTAU']
        syn.ptau = self.setting['synapse']['SYN_PTAU']
        h.setpointer(self.bcm._ref_d, 'd', syn)
        h.setpointer(self.bcm._ref_p, 'p', syn)
        return syn

    def show_synapses_PP_LTP(self):
        """Shows a box containing the model shape plot with marked recording sites and synapses."""
        # build dummy current clamps and demarcate locations
        dum_CC = [h.IClamp(self.CA1.lm_medium1(0.1)), h.IClamp(self.CA1.radTdist2(0.1))]
        for cc in dum_CC:
            cc.dur = 0
            cc.amp = 0

        self.vBoxShape = h.VBox()
        self.vBoxShape.intercept(True)
        self.shplot = h.Shape()

        for sec in self.synapses:
            for syn in self.synapses[sec]:
                if syn.type == 'AMPA':
                    self.shplot.point_mark(syn.synapse, 1, 'O', 6)

        self.shplot.point_mark(dum_CC[0], 2, 'O', 12)
        self.shplot.point_mark(dum_CC[1], 3, 'O', 12)

        self.shplot.label("Large circles: recording sites")
        self.shplot.label("Black: Synapses")
        self.shplot.exec_menu("Whole Scene")
        self.shplot.flush()
        self.vBoxShape.intercept(0)
        self.vBoxShape.map("Spatial distribution of point processes", 1200, 0, 500, 900)
        self.shplot.exec_menu("View = plot")
        self.shplot.exec_menu("Show Diam")
