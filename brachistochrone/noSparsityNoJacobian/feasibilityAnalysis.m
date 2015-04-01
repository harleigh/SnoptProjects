%{
   Feasibility analysis of the optimal solution received from Snopt to
   the Brachiscrohne problem:
 
   We use Ode45 to solve the dynamics from the Brachiscrohne problem; to
   the dynamics, we feed the optimal control received by Snopt and then
   compare the optimal solution (returned by Snopt) to the solution
   returned from Ode45.

   Dynamics:
            \dot x(t) = v(t)*sin(theta(t))
            \dot y(t) = v(t)*cos(theta(t))
            \dot v(t) = g*cos(theta(t)) - kFriction*v((t))
   where due to scaling, g=1 (so it does not appear in our dynamics).

Implementation Notes:
  * To interpolate at one point, say 0.001:
         interp1(oldMesh, theta, 0.001, 'pchip')

The following is really important:
  * We interpolate sin(theta) and cos(theta) rather than theta directly,
    this is because the optimal theta returned from the Snopt Brach problem
    is very non-smooth. Because of this, the interpolation of theta had a
    fair amount of error: Enough error to cause the feasibility to fail
    (that is, to make it look like what Snopt returned did not satisfy the
    dynamics).
  * This program tests the feasibility of the Brach Friction Problem where N
    represents N+1 sample points in the time domain [t0,tf]
%}
%clear all; close all;

%{
  We load the variables: xOptPos yOptPos vOpt tfOpt N x0 y0 v0 kFr
   * xOptPos yOptPos vOpt tfOpt are the optimal results gleaned from the
     Brach problem (located in this directory in brachFriction.m)
   * x0,y0,v0 are the initial conditions which resulted with the optimal
     results.
   * kFr is the friction coeff which resulted witht the optimal results
     above.
   * tDom is the time-domain [0,tf]
   * Note that thetaOpt is only of length N (from the first time node to the
     second to the last time node)
   * N represents N+1 sample points in tDom (that is why newMesh has an N+1 in
     it
%}
load('myVariables.mat');

%The initial point is the optimal values returned by SNOPT
xOpt0=xOptPos(1);
yOpt0=yOptPos(1);
vOpt0=vOpt(1);

%Due to the Euler Discritization Method, while the x and y dynamics are of length N, the
%theta variable is only N-1, hence we see the end-1.  See below, as we are interpolating
%based off of theta.
oldMesh = tDom(1:end-1);

%{
    Note that we are not directly interpolating the optimal theta (control) given by
    Snopt, rather we are interpolating the sin and cos of the theta-control.  The reason
    for this, is that often with this type of control, the optimal control can be very
    non-smooth due to the periodic nature of the control.  This non-smoothness can be bad
    enough that if we directly interpolated the control, we will fail feasablity due to
    error propegations.

    Note that sin(thetaOpt) contains N elements, yet tDom contains N+1 elements, hence the
    need to take off the final point in defining oldMesh above.
%}
SinThetaCtrl = @(t) interp1(oldMesh, sin(thetaOpt), t, 'pchip');
CosThetaCtrl = @(t) interp1(oldMesh, cos(thetaOpt), t, 'pchip');

%State x' = [ x y v ]
stateDynamics = @(t, x) [  x(3)*SinThetaCtrl(t); ...
                           x(3)*CosThetaCtrl(t); ...
                           CosThetaCtrl(t) - kFr*x(3)];

[t,x]=ode45(stateDynamics, [t0 tfOpt], [xOpt0,yOpt0,vOpt0]);


%%%%%%%%%%%%%%%%%%%%%%%%  Plot Code  %%%%%%%%%%%%%%%%%%%%%%%%%

%The solution returned from SNOPT is drawn with circles, and the solution
%returned from ode45 is a (blue) line.  For the solution returned by SNOPT to
%pass feasibility analysis, the solution of ode45 should match the SNOPT
%solution.
figure;
subplot(3,1,1);
    plot(t,x(:,1), tDom, xOptPos ,'o');
    title('x vs t','fontSize',14,'fontWeight','bold');
    xlabel('t');
    ylabel('x   ');
    set(get(gca,'YLabel'),'Rotation',0)
subplot(3,1,2);
    plot(t,x(:,2), tDom, yOptPos ,'o');
    title('y vs t','fontSize',14,'fontWeight','bold');
    xlabel('t');
    ylabel('y   ');
    set(get(gca,'YLabel'),'Rotation',0)
subplot(3,1,3);
    plot(t,x(:,3), tDom, vOpt ,'o');
    title('v vs t','fontSize',14,'fontWeight','bold');
    xlabel('t');
    ylabel('v   ');
    set(get(gca,'YLabel'),'Rotation',0)

%
% end feasibilityAnalysis.m
%