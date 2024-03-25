# config.py

class Config:
    def __init__(self):
        self.defaults = DefaultConfig()
        self.extreme_count_cap = self.defaults.extreme_count_cap


class DefaultConfig:
    def __init__(self):
        self.extreme_count_cap = True