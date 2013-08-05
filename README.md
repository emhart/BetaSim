Predator foraging strategies and beta diveristy
========================================================

This is the repo for the code to simulate different predator foraging preferences and estimates of beta diversity.

The initial code for the simulation was written by [Ben Bolker](http://ms.mcmaster.ca/~bolker/) and we will build upon it implementing the following forms of predations


Types of predation to implement:
-----
* Proportional: Typical generalist
  * Predator reduces all individuals by a foraging proportion which is a model parameter. Reduction fraction common across all taxa.
  
* Preferential feeding:
  * Predators eat rare things first, inverse of generalist
  * Predators prefer common things

* Predators forage with specialization:
  * Total specialization with preference for only one species per site
  * Random preference respect to abundance.

* Predator foraging different across all patches


