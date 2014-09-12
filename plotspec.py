import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt

f1 = pf.open('/Users/dm/Documents/Isotropy/art_spec/spec-0280-51612-0447.fits')
f2 = pf.open('/Users/dm/Documents/Isotropy/art_spec/spec-000001.fits')
tbdata = f2[1].data

l = tbdata['LAMBDA']
flux = tbdata['FLUX']

plt.figure()
plt.plot(l,flux)
plt.show()
