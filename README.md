# ITCZ-Breakdown-Simulation-MATLAB

Final project for MA 360 - Mathematical Modeling and Simulation  
Simulating Cyclone Formation with the Barotropic Vorticity Equation.  
Project for MA 360 Mathematical Modelling and Simulation I  

Embry-Riddle Aeronautical University  

**Created by:**  
Alexander Donato  
Samuel Rachaelson  
James Trafny  

## Background

The Barotropic Vorticity Equation (BVE) model was a method of simulating atmospheric phenomena. The program used in this experiment was built using MATLAB, with the help of meteorology professor Dr. Thomas Guinn. The goal was to simulate meso-vortices associated with a hurricane eye wall. This experiment focuses on simulating the breakdown of the Intertropical Convergence Zone (ITCZ) and subsequent formation of cyclones. The ITCZ is a line of convective thunderstorms that form near the equator and is related to the maximum incoming solar radiation from the sun. After a period of time, the ITCZ falls apart leading to the formation of tropical cyclones. This model simulates the breakdown of the ITCZ using the BVE. The goal is to explain why a particular number of storms form.

## Methodology

The previous experiment was designed to show the number of meso-vortices that formed from the breakdown of a Gaussian shape (given an initial perturbation). The breakdown occurs due to oscillating waves that manifest, feed off of each other, and phase-lock. The experiment for this project will instead be simulating the breakdown of the Intertropical Convergence Zone (ITCZ) into hurricanes. The ITCZ is an area of high activity for thunderstorms and is commonly seen in satellite imagery as a cluster of storms in a line formation as shown in Figure 1. Instead of observing the number of meso-vortices that form, we will be observing larger pools of vorticity that form. These pools of vorticity represent tropical cyclones.  

There are two primary differences between this experiment and the previous one: the shape of the feature and the scale of the feature. Instead of a Gaussian circle (resembling a Bundt cake) in order to represent a hurricane eye wall, it is a quasi-linear feature that represents the ITCZ (Sult, 2017). The scale of the feature in the previous experiment was approximately 100 km. The scale of the feature in this experiment (the ITCZ) is 2000 km.

Since the scale of this experiment is much larger, the BVE must be modified. The form of the BVE used in the previous experiment as shown below:
 
 ![Alt text](relative/path/to/img.jpg?raw=true "Title")
 
There was one term used, the Jacobian. This term represents relative vorticity advection. Vorticity is the tendency to rotate about a point in space. The model shows the evolution of areas of vorticity. There is a relationship between vorticity and the magnitude of wind speeds (the strength) of cyclones. However, it is not a direct relationship.

When dealing with a much larger feature, an extra term must be added as shown below as well as Appendix A (Page 22):
 
The second term in the equation is a derivative of a quantity, psi (known as the Stream function), with respect to x. This term represents planetary vorticity advection. It takes into consideration advection due to the rotation of the planet. Since our experiment is simulating a feature that is near the equator (9°N) and planetary vorticity advection is a function of latitude, this term only has a small effect on the simulation (Donato, 2017).

Just as with the previous experiment, an iterative numerical approach will be used to find the solution. The experiment uses a Centered-Time Centered-Space (CTCS) scheme as shown below as well as Appendix A (Page 22), Appendix B (Page 24), Appendix C (Page 25), and Appendix D (Page 26):
 
The term on the left of the equals sign represents the future iteration which requires the sum of the previous iteration and 2 multiplied by the change in t and the function evaluated at the current iteration.

This model is computationally expensive. The target length of time for the model run is two days (48 hours) since the equation governing the model is best used for the initial breakdown of the ITCZ. The model is not designed to simulate the evolution of cyclones once they form. However, longer simulations (with higher spatial and temporal resolutions) will be conducted. Animations will be created showing the evolution of the feature. The results of the experiment will be compared to a real world scenario. A satellite loop showing the ITCZ and subsequent formation of hurricanes will be used to compare to our model using the same line formation. A similar number of hurricanes should form.

## Model
The model was built in MATLAB using a pre-built engine for plotting rotation due to the BVE. The original software was heavily modified from its original form in order to achieve the testing goals. These modifications include changing the shape of the initial feature, adding an additional planetary rotation term to the BVE function, and scaling the physical grid to 2,000 km by 2000 km for the Hermite shape and 14,000 km by 40,000 km for the Earth-Scale experiments. Additional changes were made to the code base for the purposes of refactoring structural errors and memory management, as well as quality of life improvements for the output and defining grid parameters.

Appendix A (Page 20) defines the grid and time parameters. The variables jmax and imax define the grid space in the i and j directions, and thus the resolution when compared to Lx and Ly which defines the physical space. The variable dt controls the time step, decreasing this number increases the fidelity of the simulation, but also increases computation time. The variable time_hrs defines how many hours to run the simulation. 

In the Hermite version of the code, Appendix A (Page 21) defines the initial shape of the Hermite. Initially, the distance from the center line is calculated and stored into r. To create the Hermite, the grid outside of the shape was set to zero with the Hermite formula used for the sides and zeta_max for the top. The Hermite formula is shown below:
 
The Hermite equation is for each point on our grid.  
For the Earth-Scaled version (not used for validation), a simplified shape was used because of its size. For this shape, the grid was randomized and a Gaussian distribution was applied to each point centered on nine degrees latitude. The Gaussian formula is shown below:
 
The Gaussian equation is for each point on our grid, where a is the peak vorticity, b is the center of the feature (imax/2), and c is the width of the curve. Once the shape is created, the average vorticity is subtracted to start the model with an overall initial vorticity of zero.

After the shape has been created, the result is a complete zeta_grid. If it were to be plotted, Figure 7 would be the output. To get the simulation running, zeta_grid must be fed into a Jacobian matrix along with an (all zero’s) psi grid as shown in Appendix A (Page 22). The Jacobian function effectively advances our simulation by one step. Therefore the output is captured, based on the time step, and recorded to a video file. The final two “for” loops are designed to iterate the Jacobian function at a frequency of dt and record at a rate hard coded into the “if” statement. Currently it is set to 0.2 to reduce the rendering load for the project.  

The frames that are saved to the video file are then compiled into .avi format and then saved into a Results folder (which must be created in the root directory of your hard drive), programmatically named based on the temporal and spatial resolution.  It is important to make sure to change the name of any video file you wish to keep to prevent accidental overwriting.

## Validation
Our model was validated using data collected by the European Organization for the Exploitation of Meteorological Satellites (EUMETSAT) for 2017. The data uses an interpolation scheme to combine satellite imagery off the coast of South America to produce a fluid animation to fully capture an event. We used the Intertropical Convergence Zone (ITCZ) off the west coast of Panama during the month of September as our event to compare, as it produced three hurricanes during that time in close proximity to the ITCZ. Our model was set to a latitude of nine degrees north as well as 160 kilometers (100 statute miles) for the thickness of the ITCZ to resemble the event. However, the thickness of the ITCZ for the real world case varies between 80 kilometers (50 statute miles) and 320 kilometers (200 statute miles). The model ran with a spatial resolution of 128 by 128 (15.625 kilometers for grid spacing) and a temporal resolution of 10 snapshots per hour of the simulation (30 frames per second for 16 seconds for a total of 480 frames). The model was set for a 48 hour simulation. The goal was to produce a simulation with a high temporal resolution to fully capture the breakdown of the ITCZ.

## Results and Discussion 
The first 6 hours of the simulation shows no noticeable movement. However, changes are occurring in the simulation. The northern section of vorticity is moving westward and the southern section of vorticity is moving eastward. By the 12th hour, the creation of waves becomes apparent. Due to the thickness of the ITCZ, 3 defined waves form. Vorticity advection continues to amplify the waves until they phase-lock and break apart by the 24th hour. 3 cyclones form from the breakdown of the ITCZ which aligns closely with the real world case. However, the leftmost cyclone in the real world example appears to have formed before the breakdown of the ITCZ. The breakdown of the ITCZ may have supported development of the leftmost cyclone.

The last 24 hours of the simulation show the center cyclone and rightmost cyclone interacting with each other and eventually combining into one larger storm. Although cyclones can interact with each other in a real world case, direct interaction did not occur in the real world case used for this project. Changing the random seed generator for creating the disturbance of the initial vorticity field would have produced different results. However, the number of cyclones that would form from the breakdown of the ITCZ would not change since that is governed by the BVE and the thickness of the ITCZ.

## Conclusion
The BVE model was able to successfully simulate the breakdown of the ITCZ and subsequent formation of cyclones with a similar number using a similar thickness of the ITCZ compared to a real world case. This result means that the formation of 3 cyclones in the Pacific Ocean in September 2017 was determined in part by wave creation as shown in the simulation. If the thickness of the ITCZ was larger or smaller, the formation of cyclones would have changed. A larger ITCZ would results in fewer cyclones that would be potentially more powerful while a smaller ITCZ would produce more cyclones that would be potentially less powerful. If the environment was not conducive for tropical cyclone development (due to various other factors), fewer storms may have formed. The BVE model is an idealized case for the formation of cyclones. If environmental conditions are appropriate, real world formation can be linked to wave creation as was shown in Figures 1 through 3 and Figures 7 through 9.

In order to create a more realistic model, a new forcing function would need to be created that combines terms which each describe aspects of atmospheric motion (such as the pressure gradient force and frictional force). This is the approach that modern numerical models take for simulating the future state of the atmosphere. However, these models are difficult to implement and require significantly higher computational resources.
 
## References
Clay, A., Donato, A., Sult, P., & Guinn, T. (2017). BVE Model (Version 1.0) [Computer 
software]. Daytona Beach: Embry-Riddle Aeronautical University.

Green, T. & Guinn, T. (2017). relaxationcyclic.m [Computer software]. Daytona Beach: Embry-
Riddle Aeronautical University.

Geospatial Media & Communications (Producer). (n.d.). Year of Weather 2017 by EUMETSAT 
[Television broadcast]. Retrieved from https://www.youtube.com/watch?v=_2xyHIg1dA0
Hurricanes: Science and Society. (2010). Brief History of Hurricane Forecast Models. Retrieved 
from http://www.hurricanescience.org/science/forecast/models/modelshistory/

Weisstein, E. W. (n.d.). Gaussian Function. Retrieved from 
http://mathworld.wolfram.com/GaussianFunction.html

Weisstein, E. W. (n.d.). Hermite's Interpolating Polynomial. Retrieved from 
http://mathworld.wolfram.com/HermitesInterpolatingPolynomial.html

