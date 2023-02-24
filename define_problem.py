class Problem:
    def __init__(self, promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor, func=None):
        self.promo_node = promo_node
        self.max_part = max_part
        self.min_dose = min_dose
        self.max_dose = max_dose
        self.dose_interval = dose_interval
        self.inhibitor = inhibitor
        self.func = func
