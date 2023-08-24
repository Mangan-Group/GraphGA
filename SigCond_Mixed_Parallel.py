from multiprocessing import Process, Queue
import pandas as pd
from main_function_sigcond import*
import concurrent.futures
import matplotlib.pyplot as plt





def simulate(topology, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_off = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('off', topology,))[-1, -1]
    rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
    return rep_off, rep_on


def func(toplogy):
    rep_off, rep_on = simulate(toplogy)
    ON_rel = rep_on / Ref[toplogy.promo_node]['on']
    FI_rel = (rep_on / rep_off) / (Ref[toplogy.promo_node]['fi'])
    return [-ON_rel, -FI_rel]

list = []

def f(q):
    q.put([42, None, 'hello'])

# if __name__ == '__main__':
#     q = Queue()
#     p = Process(target=func, args=(q,))
#     p.start()
#     print(q.get())    # prints "[42, None, 'hello']"
#     p.join()
#
# print(list)

if __name__ == '__main__':

    with concurrent.futures.ProcessPoolExecutor() as executor:

        with open('SigCond/SigCond_combo_2.pkl', 'rb') as f:
            sigcond = pd.read_pickle(f)

        results = [executor.submit(func(topo)) for topo in sigcond]

        for f in concurrent.futures.as_completed(results):
            print(f)
