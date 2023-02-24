"""
This program is copied from pymoo github https://github.com/anyoptimization/pymoo/blob/main/pymoo/operators/survival/rank_and_crowding/classes.py.
There were errors in importing the module, therefore I had to make a copy.
"""
from rankcrowding_metrics import *
from pymoo.util.randomized_argsort import randomized_argsort
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from pymoo.core.survival import Survival #, split_by_feasibility
# from pymoo.core.population import Population


class RankAndCrowding(Survival):

    def __init__(self, nds=None, crowding_func="cd"):

        crowding_func_ = get_crowding_function(crowding_func)

        super().__init__(filter_infeasible=True)
        self.nds = nds if nds is not None else NonDominatedSorting()
        self.crowding_func = crowding_func_

    def do(self, F, n_survive=None, return_rank=False):
        # the final indices of surviving individuals
        survivors = []

        # do the non-dominated sorting until splitting front
        fronts = self.nds.do(F, n_stop_if_ranked=n_survive)

        rank_dict = dict()

        for k, front in enumerate(fronts):
            I = np.arange(len(front))
            if len(survivors) + len(I) > n_survive:
                # Define how many will be removed
                n_remove = len(survivors) + len(front) - n_survive

                # re-calculate the crowding distance of the front
                crowding_of_front = \
                    self.crowding_func.do(
                        F[front, :],
                        n_remove=n_remove
                    )

                I = randomized_argsort(crowding_of_front, order='descending', method='numpy')
                I = I[:-n_remove]
                # otherwise take the whole front unsorted
            else:
                # calculate the crowding distance of the front
                crowding_of_front = \
                    self.crowding_func.do(
                        F[front, :],
                        n_remove=0
                    )
            for j, i in enumerate(front):
                rank_dict.update({i: {"rank": k, "crowding": crowding_of_front[j]}})
            survivors.extend(front[I])
        if return_rank:
            return survivors, rank_dict
        else:
            return survivors
