The code in this repository is a simple 1D Barotropic Quasi-Geostrophic model that tracks Rossby Wave propagation.  The fundamental equation utilized in this model is the 1D potential vorticity conservation equation:

    ∂q/∂t + u ∂q/∂x + β ∂ψ/∂x = 0      

where q is the potential vorticity, u is the zonal velocity component, β is the planetary vorticity, and ψ is the streamfunction, which is related to the potential vorticity via:

    q = ∇^2 ψ 

The model utilizes a finite difference grid scheme as well as a Successive Over-Relaxation (SOR) method for time stepping (see the file: “Successive_Over-Relaxation_Method.md” for details).  

All constants are listed at the beginning of the script (“1D_Barotropic_Model_CODE.py”), and the output produces a four-panel plot that includes the following:

	Top-left: the streamfunction evolution (filled) with contour fill set to 25 levels
	Top-right: the streamfunction evolution with contours only
	Bottom-left: the streamfunction evolution (filled) with contour fill set to 50 levels
	Bottom-right: the streamfunction at a given time step (ψ vs. x)
  
An overlaid black line appears on the first three panels and is an indicator of the selected time step (the “t_index” variable).  This overlaid black line is what is displayed on the bottom-right panel that plots the streamfunction vs. the horizontal distance.  An example output is shown in the file: “Example_Output.md”

The only packages needed to run this model are numpy and matplotlib.  The goal of this model was to create a simple, easy-to-use code to examine Rossby Wave propagation.

The code here is free and available to download, use, and can be modified any way.  I do ask for an appropriate acknowledgement should it be used in any kind of project, publication, teaching material, etc.




Reference Links:

https://blogs.millersville.edu/adecaria/files/2021/11/esci342_lesson13_vorticity_equation.pdf 



