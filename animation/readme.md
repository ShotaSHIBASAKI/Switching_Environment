# Read Me
This is an example of two specie dynamics when the switching rate is $\nu=0.1$ and the txoin sensitivities of two species are $\delta=0.2$. 
We will see that time $T_{end}=200$ is enough long to converge a quasi-stationary distribution.

main.c is a code for running the dynamics. Note that in this version, the code save the dynamics at every time point certain reaction happens. The file sizes would be huge in some cases.
animation.py proviedes functions to draw gif animation. 
