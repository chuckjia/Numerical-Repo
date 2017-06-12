from Constants import *
import numpy as np

def TExactTest(x, p, t):
    return x * np.sin(t * (x + p))
