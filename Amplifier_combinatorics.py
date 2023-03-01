from Combinatorics import *
import sys
from scipy.integrate import odeint
from get_system_equations import system_equations


num_part = int(sys.argv[1])

sampling = get_combinatorics('P1', num_part, 0, 75, 0, False)
with open("Amplifier_combo_%d.pkl" % num_part, "wb") as fid:
        pickle.dump(sampling, fid)


def simulate(topology, max_time=42):
	t = np.arange(0, max_time + 1, 1)
	rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
	return rep_on

objective = [simulate(g) for g in sampling]
with open("Amplifier_objective_%d.pkl" % num_part, "wb") as fid:
	pickle.dump(objective, fid)

