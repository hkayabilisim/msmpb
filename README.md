Introduction
------------
[This Code Ocean capsule](https://codeocean.com/capsule/4293677/tree) is created to reproduce the results 
presented in the following manuscript:

H. Kaya, D. J. Hardy and R. D. Skeel, 
"Multilevel Summation for Periodic Electrostatics Using B-Splines",
Journal of Chemical Physics, submitted to Journal of Chemical Physics on 2020 December.

Please visit published Code Ocean capsule:

Mapping
-------
    Manuscript       Code Ocean Capsule
    ----------       ------------------
    Figure 1         Figure1.pdf
                     Figure1.txt (data)
    Figure 2         Figure2.pdf
                     Figure2-training, Figure2-testing-txt (data)
    Table I          Figure2.txt (fudgefactors)
    Figure 3         an illustration: not presented here 
    Figure 4         *Figure4.pdf
                     *Figure4.txt (data)
    Table III,IV     *Figure4.txt (in tabulated form)
    Table V, VI      *AppendixF.txt (in tabulated form)
    Table II         summary of Table V and VI: not presented here

    * Please read explanation about timings 


Timing
------
Please note that the timings produced in this capsule are slightly different than the ones
presented in the manuscript. Timings in the manuscript are obtained by running 
the code on a non-virtual environment with 
4 GHz Intel Core i7 CPU with 32 GB 1600 MHz DDR3 memory. Whereas the timings in the
capsule are obtained on containers running on AWS virtual machines. For this reason
they are slightly lower than the ones in the manuscript.

Code
----
* **forcefield.h, forcefield0.c, forcefield1.c**
The implementation of the algorithm presented in the manuscript.  
* **msmpb.c** 
A driver to use the above code. It is a small standalone main file.
* **pme.h, pme0.c, pme1.c**
A PME implementation to make comparions.
* **pme.c**
A driver to use PME codes. Its usage is very similar to msmpb.c.
* **msmpb_util.py**
A collection of Python wrappers to use MSM and PME.
* **AppendixF.py, Figure1.py, Figure2.py, Figure4.py**
Python codes to conduct the experiments in the manuscript. They use msmpb_util.py.
* **Makefile** 
Compiles MSM and PME source codes. 
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

Alternatively, you can export the capsule and run on your environment. 
Please consult Code Ocean documentation for further information.


