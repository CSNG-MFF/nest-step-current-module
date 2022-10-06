import nest
import numpy as np
import pandas as pd
nest.Install('stepcurrentmodule')
model = 'aeif_cond_exp'
nest.ResetKernel()
nest.SetKernelStatus({'resolution': float(0.1)})
n = nest.Create(f'{model}_sc')
t = 100.1
offset = 0
for i in range(100):
	if i > 0:
		t = 100
	I_step = 20*np.array([10., -20., 30., -40., 50.])
	t_step = np.array([0., 30., 40., 45., 50.]) + offset
	print(nest.biological_time)
	print(t_step)
	nest.SetStatus(n, {'amplitude_values':I_step, 'amplitude_times': t_step})
	nest.Simulate(t)
	offset += t

