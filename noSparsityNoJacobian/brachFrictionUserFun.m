%{
    The UserFun for the Brach Problem. This function is called by SNOPT to
    converge to an optimal solution.  In this function
        1) We pack the nonlinear constraints into F.
        2) We calculate the Jacobian of the constraints with respect to the
           decision variables (variable G).
    In this function, we discritize the dynamics of the Brach problem which were
    set as equality constraints (using Flow and Fupp in the script that calls
    this userFun).

    We discritize with respect to Euler's Method: For all i=0,...,N-1, we have
                      \dot x_{i} = (x_{i+1} - x_i)/dt
    where dt = t_{i+1}-ti == (tf-t0)/N 

    Inputs:
        decVars: Column vector whose order is determined on how xInit was
                 packed; see buildInitialGuess

    Outputs:
        F: The column vector containing the nonlinear constraints
        G: The Jacobian of the constraints with respect to the decision
           variables; G is a column vector
%}
function [ F, G ] = brachFrictionUserFun( decVars )

    global t0 kFr N D;
    
    %Unpack 'decVars' into meangful variables; memory and speed are not an issue
    %here. Hence use as many temporary variable as you want.  For performance,
    %only work with the passed-in decVars column vector.
    x = decVars(1:N+1);
    y = decVars(N+2: 2*N+2);
    v = decVars(2*N+3:3*N+3);
    theta = decVars(3*N+4: 4*N+3);
    
    %This is the objective function (final transfer time)
    tf = decVars(4*N+4);
    
    dt = (tf-t0)/N; %i.e 'h' the time step of uniform length
    
    %Here we discritize the dynamics of the Brach Problem wrt Euler's Method.
    %Note that dt is multiplied away from Euler's Differentation matrix D
    %because doing so creates a linearity in the Jacobian
    dynamicsX = D*x - dt*v(1:end-1).*sin(theta(1:end));         %colVec N by 1
    dynamicsY = D*y - dt*v(1:end-1).*cos(theta(1:end));         %colVec N by 1
    E = D; E(1:N+1:(N+1)*N)=-1 + kFr*dt;
    dynamicsV = E*v - dt*cos(theta(1:end)); %colVec N by 1
    
    %The way F is packed here determines the order of Fupp Flow and the order of
    %the Jacobian of the constraints wrt the decision variables.  F stores the
    %objective function as well as the constraints.  Note, if you set A iAfun
    %jAvar (the linear parts of the Jacobian) then do not pack those variables
    %in F--this is just a heads up; see the code example where we identify the
    %linear and nonlinear portions (and calculate the Jacobian).
    F = [ tf;
          dynamicsX;
          dynamicsY;
          dynamicsV;
          x(1);
          y(1);
          v(1);
          x(end);
          y(end) ];
     
     %The Jacobian of the constraints with respect to the decision variables; in
     %this example, we tell Snopt to do it for us.
     G=[];
end

