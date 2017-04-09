# Runge-Kutta
Runge-Kutta C program, methods (RK12 and RK24) for solving ordinary differential equations, with adaptive step size.

The program solves an ODE of the form y'(t) = f(t,y(t)) for a point M0=(t0,y(t0)) for all t in [a,b]

We propose two examples of ODE, please feel free to add or modify an ODE after defining the right conditions.

The Runge-Kutta code can be found easily on the web.
However, I decided to post my code because I had a hard time finding a RK code with adaptive step size which is very important in practice.

My code guarantees the error to be smaller than the user-defined tolerance value (it can be 10E-6, just go nuts!).
If the error is greater than the entered tolerance it means that you have entered a very small tolerance that reduced the step into the machine zero (machine epsilon). In such case, the program increases the tolerance to avoid an infinite loop.
Otherwise, if the error is MUCH smaller than the tolerance the step is augmented to optimize the program and reduce the number of points of the solution.
HINT: when the program asks for the number of steps I recommend entering (1), the program will divide the number of steps into the necessary number of points to meet the entered error tolerance. You can also enter (100000), the program reduces the number of steps as well. However, you will find that the final number of steps in both cases is a little different because the program leaves a safety marge when reducing a step.

How to use my program, create a directory with :
- main.c
- user_interface.h
- makefile
- files (create an empty folder)

Then run this directory on terminal and enter "make RK", then run the program by the command "./RK" (without the quotation marks of course!). You can simply compile and execute "main.c" but you need "-lm" to compile it for the "math.h" module.


Si vous comprenez le français, je peux vous envoyer le rapport en format PDF où j'ai bien expliqué le contrôle du pas (adaptive step-size) mathématiquement. J'ai eu la flemme de traduire le rapport en anglais...

I hope my code meets your needs.
Enjoy!
