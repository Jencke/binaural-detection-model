"""Calculate results for all experiments.

This script runs the calculation for all experiments and stores the results in
the data.h5 file. Previous data in that file is overwritten.

"""

from experiments import bernstein2014
from experiments import bernstein2020
from experiments import langford1964
from experiments import pollack1959
from experiments import rabiner1966
from experiments import robinson1963
from experiments import vandepar1999
from experiments import vanderheijden1999

# Run the experiments with default parameters
bernstein2014.calc()
pollack1959.calc()
rabiner1966.calc()
bernstein2020.calc()
langford1964.calc()
robinson1963.calc()
vandepar1999.calc()
vanderheijden1999.calc()
