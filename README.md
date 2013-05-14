10-602 Research Project
=======================

This is a research project in cosmology and machine learning. 
The goal is to learn the velocity distributions within dark matter halos using support vector regression (SVR). 
I am receiving course credit via the Machine Learning class 10-602, Independent Study: Research, at Carnegie Mellon University.

My Code
-------

*   parse_halos.py - Python parser for N-body binary output file describing halos. 
*   parse_particles.cpp - C++ parser for N-body binary output file describing individual particles.
*   constants.py - Relevant constants, including scientific quantities, simulation parameters, and named indices for use with parsed data.     		      
*   chunk_manager.py - Rough solution to allow faster manipulation of the (very large) particle data. Allows me to partition the data into .txt "chunks" such that any halo's associated particles can be accessed as quickly as the first halo's. 
*   velocity.py - Python script to generate training data and learn a parametrization of the halos' internal velocity distributions via SVR. Run "python velocity.py -h True" for more details.

Acknowledgements
----------------

My work is supervised by:

*   Shirley Ho (Assistant Professor, Physics) - shirleyh@andrew.cmu.edu
*   Jeff Schneider (Associate Research Professor, Robotics Institute / School of Computer Science) - jeff.schneider@cs.cmu.edu

I make use of dark matter N-body simulations made by:

*    Hy Trac (Assistant Professor, Physics) - hytrac@cmu.edu

I'm currently working with the following catalogs (at redshift 0):
*    halo.z=00.0000
*    subhalo.z=00.0000
*    halo_part.z=00.0000

The subhalos in subhalo.z=00.0000 were created by Michelle Ntampaka (Ph.D. student, Physics) - cmhicks@andrew.cmu.edu