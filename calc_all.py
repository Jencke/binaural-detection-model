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
print(r"--------- Calculating 1/8 Bernstein & Trahiotis 2014 ---------")
bernstein2014.calc()
print(r"--------- Calculating 2/8 Pollack & Trittioe 1959 ---------")
pollack1959.calc()
print(r"--------- Calculating 3/8 Rabiner et al. 1966 ---------")
rabiner1966.calc()
print(r"--------- Calculating 4/8 Bernstein & Trahiotis 2020 ---------")
bernstein2020.calc()
print(r"--------- Calculating 5/8 Langford & Jeffress 1964 ---------")
langford1964.calc()
print(r"--------- Calculating 6/8 Langford & Jeffress 1964 ---------")
robinson1963.calc()
print(r"--------- Calculating 7/8 van de Par & Kohlrausch 1999 ---------")
vandepar1999.calc()
print(r"--------- Calculating 8/8 van der Heijden & Trahiotis 1999 ---------")
vanderheijden1999.calc()
