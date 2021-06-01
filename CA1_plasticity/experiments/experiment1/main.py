"""
Creates CA1 pyramidal cell model, sets parameters, runs simulation in NEURON and shows figures.

Title: main.py
Author: Matus Tomko
Mail: matus.tomko __at__ fmph.uniba.sk
"""
from CA1_plasticity.io.figures_shower import FiguresShower
from CA1_plasticity.io.io_helper import IOHelper
from CA1_plasticity.model.CA1_pyramidal_cell import CA1PyramidalCell
from CA1_plasticity.model.stimulation_protocols import StimulationProtocol


def prepare_simulation():
    CA1_cell = CA1PyramidalCell(hoc_model='../../model/hoc_models/pc_strong_bAP_updated.hoc',
                                path_mods='../../model/mods/', setting=setting)
    CA1_cell.insert_AP_counter()
    CA1_cell.insert_BCM()
    CA1_cell.insert_json_synapses(synapses=ioh.load_synapses())
    CA1_cell.connect_spontaneous_activity(lm_syn_count=None, ori_syn_count=None, rad_syn_count=None,
                                          random_weights=False)

    sp = StimulationProtocol(setting=setting)
    sp.set_Dong_sequential_stimulation(CA1_cell.synapses)

    run_simulation()

    ioh.save_recordings(synapses=CA1_cell.synapses, tw_vec=CA1_cell.t_rs_vec, v_soma_vec=CA1_cell.v_vec,
                        t_vec=CA1_cell.t_vec, p_vec=CA1_cell.p_vec, d_vec=CA1_cell.d_vec, ta_vec=CA1_cell.t_rs_vec,
                        alpha_scount_vec=CA1_cell.alpha_scout_vec, dend_vecs=CA1_cell.dend_vecs,
                        apc_vec=CA1_cell.apc_vec, cai_vecs=CA1_cell.cai_vecs, cal_ica_vecs=CA1_cell.cal2_ica_vecs,
                        nmda_ica_vecs=CA1_cell.nmda_ica_vecs, ina_vecs=CA1_cell.ina_vecs)
    ioh.save_setting(setting=setting)
    ioh.save_synapses(synapses=CA1_cell.synapses)


def run_simulation():
    from neuron import h, gui
    h.tstop = setting['simulation']['TSTOP']
    h.dt = setting['simulation']['DT']
    h.v_init = -65
    h.celsius = 35
    h.finitialize(-65)
    h.fcurrent()
    h.cvode_active(0)
    h.run()
    del h, gui


def figures():
    fs = FiguresShower(setting=setting, save_figures=False, path_saving='figures/example/',
                       path_recordings='recordings/example/')
    fs.show_somatic_voltage()
    fs.show_input_spikes()
    fs.show_SCH_COM_weights()
    fs.show_d_p_amplitude()
    fs.show_alpha_scount()


if __name__ == '__main__':
    ioh = IOHelper(path_saving='recordings/example/', path_settings='settings/')
    try:
        setting = ioh.load_setting()
    except Exception as e:
        print(e)
        print('Loading setting failed. Program terminated.')
        prepare_simulation()
        figures()
