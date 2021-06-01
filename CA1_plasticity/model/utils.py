"""
Title: utils.py
Author: Matus Tomko
Mail: matus.tomko __at__ fmph.uniba.sk
"""


class RecordingVector:
    """
    A class used to represent a vector for recording voltage or current

    ...

    Attributes
    ----------
    section : str
        the section name
    segment_x : float
        the recording site on the section
    vector : neuron.hoc.HocObject
        the vector for recording voltage or current
    """

    def __init__(self, section, segment_x, vec):
        """
        Parameters
        ----------
        section : str
            the section name
        segment_x : float
            the recording site on the section
        vec : neuron.hoc.HocObject
            the vector for recording voltage or current
        """

        self.section = section
        self.segment_x = segment_x
        self.vector = vec


class Synapse:
    """
    A class used to represent a Synapse

    ...

    Attributes
    ----------
    synapse : neuron.hoc.HocObject
        the HocObject of the synapse
    section : str
            the section name
    segment_x : float
        the position of the synapse on the section
    distance : float
        the distance of the synapse from the soma
    init_weight : float
        the initial weight of the synapse
    weight_vec : neuron.hoc.HocObject
        the vector for recording the synaptic weight over time
    input_spikes_vec : neuron.hoc.HocObject
        the vector for recording the times of the spikes arriving to the synapse
    stimulated : bool, optional
        the attribute determining the application of the stimulation protocol to the synapse (default False)
    receptor : str
        the receptor of the synapse (AMPA, NMDA, or None)
    pathway : str, optional
        the hippocampal pathway representing by the synapse (COM, SCH, LM, or None)
    """

    def __init__(self, synapse, section, segment_x, distance, init_weight, weight_vec, input_spikes_vec, receptor):
        """
        Parameters
        ----------
        synapse : neuron.hoc.HocObject
            the HocObject of the synapse
        section : str
            the section name
        segment_x : float
            the position of the synapse on the section
        distance : float
            the distance of the synapse from the soma
        init_weight : float
            the initial weight of the synapse
        weight_vec : neuron.hoc.HocObject
            the vector for recording the synaptic weight over time
        input_spikes_vec : neuron.hoc.HocObject
            the vector for recording the times of the spikes arriving to the synapse
        receptor : str
            the receptor of the synapse (AMPA, NMDA, or None)
        """

        self.synapse = synapse
        self.section = section
        self.segment_x = segment_x
        self.distance = distance
        self.init_weight = init_weight
        self.weight_vec = weight_vec
        self.input_spikes = input_spikes_vec
        self.stimulated = False
        self.receptor = receptor
        self.pathway = None
