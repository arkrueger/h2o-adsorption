# h2o-adsorption
Numerical solution for modeling H2O adsorption on a silica column. Fitting process parameters using finite difference and nonlinear regression.

I originally did this project in MATLAB but rewrote it in Python for my sanity.

The adsorption of water vapor onto silica can be modeled by a system of partial differential equations, whose critical parameter is the mass transfer coefficient, kappa.

Here is a brief overview of the analysis.

Step 0: Define the problem. What is the mass transfer coefficient (kappa) for water adsorption on silica in a packed column?

Step 1: Pack a column with color-changing silica beads (orange -> white).

Step 3: Flow humidified air through the silica column. At the same time, record video of the column as it changes color.

Step 4: Stop recording once the column has completed the color change.

Step 5: Analyze the images and extract the blue values. That's the B in RGB, and what it really tells me here is whether the column is orange or white.

Step 6: For each image, take slices of the column perpendicular to the direction of gas flow, average the B value, and normalize to the brightest portion of the column. In other words, how far has the water saturation progressed up the column?

Step 7: Time for math. Construct a finite difference model from a system of partial differential equations, which allows us to skip the analytical solution and solve by curve fitting. In nonlinear regression, a curve is constructed from a "guess", compared against experimental data, and scored based on goodness of fit. This continues iteratively with new guesses until the algorithm converges on the best fit. In this case, the guess is kappa. The guessed curve is generated with the aforementioned finite difference solution, and its goodness of fit is calculated by the sum of squared differences. This process is repeated with different guesses of kappa until the sum of squared differences is sufficiently minimized, i.e. the guessed kappa value produces a theoretical curve that is very close in profile to the experimentally generated data.

Step 8: Reflect on the pain caused by Step 7. Sticks and stones _may_ break my bones but MATLAB will _always_ hurt me.

Step 9: Rewrite the project in Python.

![image](https://user-images.githubusercontent.com/54046534/115321378-87971c80-a138-11eb-819c-32dc539610d7.png)
![image](https://user-images.githubusercontent.com/54046534/115322474-fecdb000-a13a-11eb-8b84-97cda945d79a.png)



And here's a snippet of the original MATLAB code before I rewrote it in Python. 
![image](https://user-images.githubusercontent.com/54046534/115321300-633b4000-a138-11eb-90e1-310993e43c87.png)

MATLAB, never again.
