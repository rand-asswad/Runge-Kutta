# Runge-Kutta
Runge-Kutta C program, methods (RK12 and RK24) for solving ordinary differential equations, with adaptive step size.

The program solves an ODE of the form y'(t) = f(t,y(t)) for a point M0=(t0,y(t0)) for all t in [a,b]

We propose two examples of ODE, please feel free to add or modify an ODE which you check its right conditions.

The Runge-Kutta code can be found easily on the web.
However, I decided to post my code because I had a hard time finding a RK code with adaptive step size which is very important in practice.

My code guarantees the error to be smaller than the user-defined tolerance value (it can be 10E-6, just go nuts!).
If the error is greater than the entered tolerance it means that you have entered a very small tolerance that reduced the step into the machine zero (machine epsilon),
in such case, the program increases the tolerance to avoid an infinite loop.
Otherwise, when the error is MUCH smaller thant the tolerance the step is augmented to optimize the program and reduce the number of points of the solution.

How to use my program, create a directory with :
- main.c
- user_interface.h
- files (create an empty folder)

I hope my code meets your needs.
Enjoy!
