# SnoptProjects
Optimal Control Problems implemented in Snopt within MatLab 

Contained within this directory are various Optimal Control
problems numerically solved using the MatLab interface to the
nonlinear solver SNOPT.

These projects will run on the student/trial versions of SNOPT,
which can be obtained at http://www.scicomp.ucsd.edu/~peg/Software.html

Projects contained within this directory:
  * Brachistochrone: Here Snopt is used to numerically solve the
    brachistochrone problem under the influence of a resistive force
    proportional to the velocity of the particle.  Various implementations
    are given, namely showcasing the explicitly computation of the Jacobian
    (of the decision variables with respect to the constraints) to increase
    convergence to the optimal solution, along with the computation of the
    sparsity pattern of the Jacobian.  Similarly, a version where Snopt
    calculates the Jacobian is given; comparing the two shows the order-of-
    magnitude in run-time.  Included in each is a feasibility analysis to
    the solution Snopt returns.  Feasibility Analysis is a critical component
    to verifying and validating that the Numerical Solver (in this case SNOPT)
    indeed has found an optimal solution.
