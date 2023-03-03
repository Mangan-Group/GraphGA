from amplifier_problem import Amplifier
from define_circuit import Topo

amp_testcase = Amplifier(
    'P1',
    1,
    75, 
    75, 
    5, 
    False, 
    {}, 
    1, 
    pop=True,
    num_processes=8
    )


manual_topo = Topo([('P1', 'Z6'), ('Z6', 'Rep')], {'Z6': 75}, 'P1')


if __name__ == '__main__':

    manual_ON_rel = amp_testcase.func(manual_topo)
    print(manual_ON_rel)
