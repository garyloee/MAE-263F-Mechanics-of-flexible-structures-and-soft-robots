Load in GradEs, GradEb, HessEs, HessEb to allow the codes to run properly.

Problem 1
run Problem1.m to excute implicit method for question 1 and 2.

To set the radius of the three spheres to the same value, change R2 on row 14.

Implicit model's timestep can also be changed by typing in the new dt value on row 7.

The explicit model will run automatically by the finction problem1Explicit.m after the implicit is finished where we can see the shape, position, velocity of the structure.

To compare different time step size, change the dt value in problem1Explicit.m to see the difference.




Problem 2
Run Problem2.m to see the solution for question 1 and 2.

To compare the effect of spatial discretization, the function problem2Spat will run after the first code. Where we will see the shape of structure getting finer nad the terminal velocity displayed in command window.

The Terminal velocity vs spacial discretization graph will also pop up after the simulation finished.

Similarly, the problem2Temp will execute and create graph to demonstrate time step size's effect on terminal velocity.

To avoid interference and save time, the function calling on row 259 and 260 of Problem2.m can be commented.



Problem 3
Execute Problem3.m to plot the infomation for question 1.

The function problem3Euler will run afterward and plot the comparison between euler beam function and our simulation.






