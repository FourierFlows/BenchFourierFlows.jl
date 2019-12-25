import numpy as np
from numpy import pi
import pyqg
import time
from pyqg.diagnostic_tools import spec_var

# the model object
m = pyqg.BTModel(L=2.*pi, nx=512, tmax=5.,
        beta=0., H=1., rek=0., rd=None, dt=0.01,
        taveint=5., ntd=4)

# McWilliams 84 IC condition
fk = m.wv != 0
ckappa = np.zeros_like(m.wv2)
ckappa[fk] = np.sqrt( m.wv2[fk]*(1. + (m.wv2[fk]/36.)**2) )**-1

nhx, nhy = m.wv2.shape

Pi_hat = np.random.randn(nhx, nhy)*ckappa +1j*np.random.randn(nhx,nhy)*ckappa

Pi = m.ifft( Pi_hat[np.newaxis] )
Pi = Pi - Pi.mean()
Pi_hat = m.fft( Pi )
KEaux = spec_var( m, m.filtr*m.wv*Pi_hat )

pih = ( Pi_hat * np.sqrt(0.5/KEaux) )
qih = -m.wv2*pih
qi = m.ifft(qih)
m.set_q(qi)



t = time.time()
# run the model
m.run()
# do stuff
elapsed = time.time() - t
print(elapsed)