"""Plot results from all experiments.

This script loads simulate results from the data.h5 file and plots model
results vs expermental results it in a 3x3 raster

"""

from experiments import bernstein2014
from experiments import bernstein2020
from experiments import langford1964
from experiments import pollack1959
from experiments import rabiner1966
from experiments import robinson1963
from experiments import vandepar1999
from experiments import vanderheijden1999

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

plt.ioff()
# Create a 3x3 raster to contain the plots
fig = plt.figure(figsize=(15, 8))
gs = GridSpec(3, 3, bottom=0.07, top=0.95, left=0.07, right=0.98,
              hspace=0.5, wspace=0.3)
sub_gs = [gs[i, j] for i in range(3) for j in range(3)]

bernstein2014.plot(sub_gs[0])
pollack1959.plot(sub_gs[1])
rabiner1966.plot(sub_gs[2])
bernstein2020.plot(sub_gs[3])
langford1964.plot(sub_gs[4])
robinson1963.plot(sub_gs[5])
vandepar1999.plot(sub_gs[6])
vanderheijden1999.plot(sub_gs[7])

plt.show()
