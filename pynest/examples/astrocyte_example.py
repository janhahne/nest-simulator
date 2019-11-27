import numpy as np
from matplotlib import pyplot as plt
import nest

nrn_model = 'aeif_cond_alpha_astro'
astro_model = 'astrocyte'
conn_nrn_astro = 'tsodyks_synapse'
conn_nrn_nrn = 'tsodyks_synapse'
conn_astro_nrn = 'sic_connection'

send_nrn = nest.Create(nrn_model, 1)
rec_nrn = nest.Create(nrn_model, 1)
astro = nest.Create(astro_model, 1)
mm = nest.Create('multimeter', params = {'record_from': ['V_m'], 'withtime': True})
curr_gen = nest.Create('dc_generator', params = {'start':1.0, 'amplitude': 800.0})

nest.Connect(send_nrn, rec_nrn, syn_spec = {'model': conn_nrn_nrn})
nest.Connect(send_nrn, astro, syn_spec = {'model': conn_nrn_astro})
nest.Connect(curr_gen, send_nrn)
#nest.Connect(astro, rec_nrn, syn_spec = {'model': conn_astro_nrn})
nest.Connect(mm, send_nrn)

nest.Simulate(150.0)

mm_data = nest.GetStatus(mm)[0]['events']
times = mm_data['times']
voltages = mm_data['V_m']

plt.figure()
plt.plot(times, voltages)
plt.show()
