from amplifier_problem import Amplifier
from define_circuit import Topo

amp_testcase = Amplifier(
    'P1',
    [75, 75, 5],
    2, 
    False, 
    False,
    {}, 
    1, 
    True,
    8
    )


manual_topo = Topo(
    # [('P1', 'Z1'),('Z1', 'Rep')],

    # [('P1', 'Z2'), ('Z6', 'Rep'), ('Z6', 'Z6'), ('Z2', 'Z6')],
    [('P1', 'Z2'), ('Z2', 'Z6'), ('Z6', 'Z6'), ('Z6', 'Rep')],
    {'Z6': 75, 'Z2': 75}, 'P1')


manual_ON_rel = amp_testcase.func(manual_topo)
print(manual_ON_rel)
# print(manual_topo.in_dict.keys())
# manual_topo.plot_graph()
