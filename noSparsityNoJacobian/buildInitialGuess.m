%{
    Name: buildInitialGuess

    This function returns an initial guess (of a straight line) that starts
    at (x0,y0) and ends at (xf,yf) and has constant velocity, and a final time
    guess of 1.

    Inputs:
        * x0, y0, v0, xf, yf Initial and final point constraints for the Brach
          Curve
        * N is N+1 uniformly distributed points on [t0,tf]
    
    Output: Initial guess as a column vector of length 4*N+4, packed as
       [x0,x1,...,xN,  y0,y1,...,yN,  v0,v1,...,vN,  theta0,...,thetaN-1,  tf]'
%}
function [ xInit ] = buildInitialGuess(x0, y0, v0, xf, yf, N)

    %The following code makes the initial guess a straight line with
    %constant velocity (of v0) and constant theta (of 0 degrees)
    
    %We find N+1 equally spaced sample points for x y v and  only N for theta;
    %this is due to the discritization scheme.  We have that x y
    %v are evaluated at equally spaced t0, t1, t2, ... , tN (Hence N+1
    %sample points).  Notation: x(ti) is denoted xi for i=0,1,...,N  (as is for the others).
    x = linspace(x0, xf, N+1); %dim 1xN+1 [x0,x1,...,xN]
    y = linspace(y0, yf, N+1); %dim 1xN+1 [y0,y1,...,yN]
    v = linspace(v0, v0, N+1); %dim 1*N+1 [v0,v1,...,vN]
    
    %As we are applying Euler's Discritization, we only need N theta points.  That is, in
    %the calculation of the dynamics, thetaN is never used, so we never define it.  If we
    %defined a decsision variable that we never use, we make our problem not well-defined.
    theta  = linspace(0*pi/4, 0, N); %dim 1xN [theta0,...,thetaN-1]
    
    %A guess on when the particle will land at (xf,yf) when at time zero the
    %particle mass starts at (x0, y0)
    tf = 1;
    
    %pack the initial guess as a column vector for SNOPT
    xInit = [ x y v theta tf ]';
    
end

