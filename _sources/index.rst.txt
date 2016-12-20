.. growthsim documentation master file, created by
   sphinx-quickstart on Wed Dec 14 18:37:02 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to GrowthSim's documentation!
=====================================

--------
Overview
--------

GrowthSim is a Python package designed to test the metabolic load
of a synthetic genetic circuit in a cell by simulating the growth
of multiple cell strains in the same environment.  The model framework
constrains protein production with the availability of nutrients, which
are made available in a resource pool shared by all cells.  This allows
different strains to outcompete one another over time.  The package
provides tools to build cell strains with customized proteins and interaction
pathways, merge a set of strains into a full population, and simulate
the population's growth over time with user-defined initial conditions.

-------------------
Classes and Methods
-------------------

.. class:: strain(self)

   The strain class extends the sbmlModel class from the simplesbml module, with
   additional functions for adding protein pathways and creating deep copies of itself.
   Upon initialization, it automatically loads a model of nutrient transport, metabolism,
   transcription, ribosome-RNA binding, translation, and degradation.  Only one cell is
   tracked, and the external nutrient concentration is held constant, but this is changed
   when the strain is added to a population.

   .. function:: addProtein(self, protein_id, n=300, txn_rate='we*a/(thetax + a)', modifiers=['a'], deg_rate='0')

      Add a protein with a given ID to the strain.  *n* is the protein length.
      *txn_rate* is an expression describing the protein's transcription rate. Add
      parameters to the strain using addParameter() before calling addProtein().
      Also remember to include a rate-limiting term for *a* in Michaelis-Menten form.
      *modifiers* is a list of species used in the transcription rate law.
      *deg_rate* is an expression for protein degradation rate, which defaults to 0.

   .. function:: clone(self)

          Create a deep copy of the current strain and return it.

.. class:: population(self, tag_dict)

   The population class also extends simplesbml.sbmlModel.  It has an additional attribute
   strains, which is a dictionary of strain objects paired with "ID tag" strings as keys.
   Upon initialization, the strains are processed to set appropriate initial conditions and
   be tracked as a population of cells, with the population size as a state variable.  External
   nutrient concentration becomes a state variable with influx and efflux, simulating a chemostat.

   .. function:: addStrain(self, modn, tag)

   Process and add a new strain to the population.  The tag argument is a string that is appended
   to the name of each state variable in the strain to distinguish it from other strains.

   .. function:: simulate(self, tstart, tend, nsteps)

   Simulate the population from time *tstart* to time *tend* with a total of *nsteps* timesteps.
   Returns a list of the population size for each strain at each timestep.

--------
Examples
--------

Here is an example of a population with two variants of a strain with a repressilator
(without a reporter protein) and a "wild-type" strain with no additional proteins::

	 import numpy as np
	 import matplotlib.pyplot as plt
	 import growthsim

	 model = growthsim.strain();
	 model.addParameter("wg", 100);
	 model.addParameter("Kg", 100);
	 model.addParameter("hg", 2);
	 model.addParameter("dg", np.log(2)/4);
	 model.addProtein("g1", 300, "wg * a / (thetax + a) / (1 + pow(g3/Kg, hg))", \
         ["a", "g3"], "dg");
	 model.addProtein("g2", 300, "wg * a / (thetax + a) / (1 + pow(g1/Kg, hg))", \
      	 ["a", "g1"], "dg");
	 model.addProtein("g3", 300, "wg * a / (thetax + a) / (1 + pow(g2/Kg, hg))", \
         ["a", "g2"], "dg");

	 model2 = growthsim.strain();

	 model3 = model.clone();
	 model3.getParameter("wg").setValue(1);

	 tag_dict = {'anc':model2.clone(), 'rep':model.clone(), 'knock':model3.clone()};
	 model_combo = growthsim.population(tag_dict);

	 result = model_combo.simulate(0, 10000, 1001);

	 labels = sorted(model_combo.strains.keys());
	 for i in range(len(labels)):
    	     plt.plot(result[:,0], result[:,i+1], label="N_" + labels[i]);
	 plt.legend();

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
