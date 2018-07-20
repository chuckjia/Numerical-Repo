#!/usr/bin/env python2.7

import numpy as np
from atm import *
from mlplot import *


[t,Nx,Np,X,Z,T,q,u,w]=read_txt('run_physic3/DATA/0.data')


plot_contour(1,T,q,t,0,X,Z,X,Z,[],[],'T_q','T_q',1,0)
