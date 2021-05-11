Introduction
------------
This GitHub [repository](https://github.com/hkayabilisim/msmpb) contains up-date codes to implement
the algorithms published in:

H. Kaya, D. J. Hardy and R. D. Skeel, 
"Multilevel Summation for Periodic Electrostatics Using B-Splines",
The Journal of Chemical Physics 154(14):144105

An older version of the code we used at the publication time is also provided 
as [a Code Ocean capsule](https://codeocean.com/capsule/4293677/tree) for those
reader who want to recreate the figures and tables in the manuscript.

For more information please consult [Wiki](https://github.com/hkayabilisim/msmpb/wiki) page.

Code
----
* **forcefield.h, forcefield0.c, forcefield1.c**
The implementation of the proposed algorithm.
* **example.c** and **h20.dat**
A simple driver and sample data for demonstration purposes.
* **msmpb.c** 
Another driver with more features. 
* **pme.h, pme0.c, pme1.c**
A PME implementation to make comparions.
* **pme.c**
A driver to use PME implementation. Its usage is very similar to msmpb.c.
* **journal** folder contains the scripts to produce the experiments in the journal paper.
* **msmpb_util.py**
A collection of Python wrappers to use MSM and PME.
* **Makefile** 
Compiles MSM and PME source codes and the drivers.
* **run** 
An entrypoint for the Code Ocean Capsule. 

Data
----
There are three different files in the data folder: *.ini, *.pot and *.acc files.
ini files contain the charges, positions and geometry. pot and acc files
contain potential energy and accelerations in high precision. spce files correspond 
to the NIST SPC/E water benchmarks. CsCl and NaCl contain the crystal structures. 
Lastly Box44 is another water benchmark.

Reproducing the Results
-----------------------
A published Code Ocean Capsule is already verified to be reproducible.
If you also want to run the capsule one more time and get the results,
you can use the tools provided in the Code Ocean's Web Site.  
Please consult Code Ocean documentation for further information.

Alternatively, you can run the scripts under the journal folder
on your local environment.

