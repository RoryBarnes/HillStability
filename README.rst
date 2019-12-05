Hill Stability of a Three-Body System
=========================

This code calculates the relative proximity of a three-body system to "Hill stability." A Hill stable system is one for which the ordering of the
bodies remains constant, i.e. the most distant body may escape to infinity and
the system would still be considered "Hill stable." Hill stability may be calculated for any three-body system, but this code is optimized for planetary systems and was first used in `Barnes, R. &
Greenberg, R. (2006) <https://ui.adsabs.harvard.edu/abs/2006ApJ...647L.163B/abstract>`_.

To compile:

.. code-block:: bash

  gcc -o hillstab hill_stab.c -lm

To execute the example:

.. code-block:: bash

  hillstab jupsat.in

This will calculate 2 types of Hill stability: Exact and Approximate. Both
numbers represent relative proximity to the boundary, with unity on the
boundary, values < 1 are Hill unstable, and > 1 are Hill stable. Exact is the
value of beta/beta_{crit} from BG06, which is computed by calculating the energy
and angular momentum including the primary. "Approx" is the value of delta/delta_{crit} from Barnes & Greenberg (2006), which is computed assuming the central mass
dominates, see Gladman (1993). As this option assumes the central body's center
is very close to the system's center-of-mass, Approx uses the input elements.   

The input file should have the following format:

.. code-block:: bash

    CentralMass
    Mass SemiMajorAxis Eccentricity Inclination ArgPeri LongAscNode MeanAnomaly
    Mass SemiMajorAxis Eccentricity Inclination ArgPeri LongAscNode MeanAnomaly
    CoordinateSystem

where the CentralMass units are in solar masses, and the following two lines
represent the two orbiters. Their units are Jupiter masses, AU, and degrees. The final line must state either "bodycentric" or "barycentric" to indicate the coordinate system of the orbital elements. There are no command line options.

If you use this code, please cite `Barnes &
Greenberg (2006) <https://ui.adsabs.harvard.edu/abs/2006ApJ...647L.163B/abstract>`_.
