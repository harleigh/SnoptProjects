%{
    What we are doing:
        A point mass is allowed to slide over riged two dimensional surface with
        (or without) a friction force acting oppisite to the velocity vector
        (assuming the friction law is given as \mu*N where \mu is the coeff of
        sliding friction) along with the influance of gravity: What is the shape
        of the surface that makes the point mass to move from one given point to
        another given point (initial and final point) in Minimum Time

    Formulation of the Brach Problem:
        Minimize   tf
        subject to the dynamics
             \dot x - v*sin(theta) = 0
             \dot y - v*cos(theta) = 0
             \dot v - g*cos(theta) - kFriction*v = 0
        subject to the endpoint constraints
             t0  = 0
           x(t0) = x0   i.e.: the initial position in x
           y(t0) = y0   i.e.: the initial position in y
           v(t0) = v0   i.e.: the initial velocity of the point mass
           x(tf) = xf   i.e.: the final position of x
           y(tf) = yf   i.e.: the final position of y

    Notes:
       * We take g to be 1 by scaling v.
         Since g is measured in m/s^2 we can apply different units
         and hence scale the problem to make g=1. For instance, 
         take \hat v = (1/g)v, then our velocity dynamics are
         \dot\hat v - cos(theta) - kFriction*cos(theta)=0
         Hence throughout the problem, we just take g=1

    In this version, N represents N+1 sample points on [t0,tf]
    and hence N intervals.
        eg: N=4 gives t0--t1--t2--t3--t4 with four intervals

    In this version, we supply the sparsity pattern of the Jacobian (of the
    constraints wrt the decision variables) as well as explicitly calculate
    every entry of the Jacobian.
  
    Note that in this version we do not supply the linear portion of the
    Jacobian; that is, we treat everything as nonlinear.

%}
clear all; close all;

%{
    Global variables nessisary to persist variables into the UserFun that Snopt
    uses to converge to the optimal solution.
        t0: Initial Time (it's always zero)
        kFr:     Kenetic friction
        N:       the N+1 sample Nodes on [t0,tf] where the first node is t0 and
                 the (N+1)th node is tf
        D:       Euler's Differentation Matrix
        index4G: This is a technique; Snopt wants the Jacobian to be a column
                 vector (which is a real pain) so index4G is used to turn the
                 computed jacobian (by me) into a column vector. (the idea is
                 that when we build the explicit Jacobian, we do so as a full
                 sized matrix (which simplifies the indexing) and when we are
                 done, we pull out the nonzero elements from the full-sized
                 matrix in the same order as the sparsity (of the Jacobi) was
                 specified. What a pain! :)
%}
global t0 kFr N D index4G;

%Set initial and final Conditions to the Brach problem.
x0 =  0; y0 = 0;   % initial position. 
xf = 5; yf = 2;   % final position.   
v0 =  sqrt(2*y0);  % initial velocity.
%--end setting initial and final Conditions

%Set simulation parameters.
t0 = 0;      % starting time
kFr = 0.0;   % kenetic Friction constant (0.6 is a good one)
N = 50;      % Num of intervals on [t0,tf], hence N+1 sample points on [t0,tf]

%This is the upper bound for the final time to the Brach problem.  Note that this value is
%a guess.  If the guess is too small, snopt might complain that the problem is infeasible
%(DIDO does), so this upper bound is a thing to always keep in mind when looking at your
%results.  If your states/costates are stopping at the upper bound for the final time, yet
%are not feasible, try increasing this bound and re-run the simulation.
FinalTimeUpperBound = 20;

%Build the Euler Differentation Matrix.  Note: the step h is multiplied through
%to the dymanics to case part of the Jacobian to be linear :D --see my notes.
D = zeros(N,N+1);
D(1:N+1:(N+1)*N)=-1;
D(N+1:N+1:(N+1)*N)=1;
%--end setting the simulation parameters


%pointer to the snopt userfun; called in the cmex to snopt.
usrfun = 'brachFrictionUserFun';

%Note, for snopt, the initial guess needs to be as a column vector.
xInit = buildInitialGuess(x0, y0, v0, xf, yf, N);

numDecVars = length(xInit);

%{
   Set the lower and upper bounds on the decision variables. SNOPT wants xLow
   and xUpp as column vectors.  Note that the order in which xInit was pack
   determines the ordering of the decision variables (and as final time was
   packed last, the last entry for xLow and xUpp are the lower and upper bounds
   for the final time).
%}
xLow = -(1.0e+06)*ones(numDecVars,1);
xUpp =  (1.0e+06)*ones(numDecVars,1);
xLow(end) = t0 ;    %lower bound for final-time.
xUpp(end) = FinalTimeUpperBound;    %upper bound for final-time.

%{
   Set the lower and upper bounds for the cost function and the constraints
   Note: the order of Flow and Fupp is determined by the packing of the column
   vector F, which is done in the userFun.  F is packed as:
        cost function--inequality constraint
        discritiezed dynamics for x position--equality constraint
        discritiezed dynamics for y position--equality constraint
        discritiezed dynamics for v position--equality constraint
        initial x of current solution
        initial y of current solution
        initial v of current solution
        final x of current solution
        final y of current solution
%}
Flow = [ t0;            %lowB bound of the cost (final time)
         zeros(3*N,1);  %lowB of the dynamics of x y and v (each size N)
         x0;            %equality constraint for initial x position of particle
         y0;            %equality constraint for initial y position of particle
         v0;            %equality constraint for initial v of particle 
         xf;            %equality constraint for final x position of particle
         yf ];          %equality constraint for final y position of particle

Fupp = [ FinalTimeUpperBound;
         zeros(3*N,1);
         x0;
         y0;
         v0;
         xf;
         yf ];
     
%This little 4-line block of setting variables for Snopt are not used very
%often; we set them to the defaults; ObjRow for example defines which row the
%objective function is located within the constraints col-vec ('F' in the
%userFun)
numConstraints  = length(Flow);
xMul   = zeros(numDecVars,1); xState = zeros(numDecVars,1);
Fmul   = zeros(numConstraints,1); Fstate = zeros(numConstraints,1);
ObjAdd = 0; ObjRow = 1;

%Setting the linear and nonlinear sparsity patterns--this version does not
%supply any of the linear sparsity patterns; these sparsity variables tell Snopt
%which elements of the Jacobian need to be calculated, and which entries can be
%ignored
A  = [];  iAfun = [];  jAvar = [];                  %linear
jacobianSparsity  = findSparsityPattern( D, N );    %nonlinear
[iGfun,jGvar] = find(jacobianSparsity); %row and col indicies of nonzero elmnts
index4G = find(jacobianSparsity);       %logical indexing of nonzero elmnts


%Set the Optimal Parameters for SNOPT. See chapter 7, pg 62 'Optimal Parameters'
%Note we first set 'Defaults' to start SNOPT on a clean-slate; very important!
snset('Defaults');              %Needed to flush Snopt; without this Snopt can Snot
snseti('Derivative option',1);  %This time we have all of the derivatives hence '1'
snseti('Verify level', 3);

%Sumary and Solution files; see chapter 8 of SNOPT guide (section 8.8, 8.9)
snprint('resultGradySNOPT.txt');
snsummary('resultGradySNOPT.sum');

%Call snopt
solveopt = 1;   %not sure what this does other than say 'Hey optimize!'
tic             %time the execution of snopt
[xOpt,F,xMul,Fmul,INFO] = snoptcmex( solveopt, ...
				                     xInit,xLow,xUpp,xMul,xState, ...
				                     Flow,Fupp,Fmul,Fstate,        ...
				                     ObjAdd,ObjRow,A,iAfun,jAvar,  ...
				                     iGfun,jGvar,usrfun );
runTime=toc;
%%%%%%%%%%%%%%%%%%%%%% Plot Code %%%%%%%%%%%%%%%%%%%%%%

snsummary off; % close the summary .sum file; empties the print buffer
snprint   off; % Closes the print-file .out; empties the print buffer

%Unpack the optimal solution from the column vector Snopt returned
xOptPos     = xOpt(1:N+1);          %optimal x        length N+1
yOptPos     = xOpt(N+2: 2*N+2);     %optimal y        length N+1
vOpt        = xOpt(2*N+3:3*N+3);    %optimal v        length N+1
thetaOpt    = xOpt(3*N+4: 4*N+3);   %optimal control  length N
tfOpt       = xOpt(4*N+4);          %final time       length 1

%unpack the initial positions (just to plot the initial guess along with the
%optimal guess
xInitPos = xInit(1:N+1);
yInitPos = xInit(N+2:2*N+2);

%Just some nicely formated output to the console on the Snopt Run
disp(strcat('Execution time: ', num2str(runTime)));
disp(strcat('SNOPT exited with INFO==', num2str(INFO)));    %see pg 19 for INFO
disp(strcat('Final Time ', num2str(tfOpt)));
disp(strcat('Initial Conditions: (x0,y0, v0)==(', num2str(x0),',',num2str(y0),',',num2str(v0),')'));
disp(strcat('Final Conditions: (xf,yf)==(', num2str(xf),',',num2str(yf),')'));
disp(strcat('N==',num2str(N), '  (hence,', num2str(N+1),' sample points on [t0,tf])'));
disp(strcat('Kenetic Friction kFr==', num2str(kFr)));

%Plot the optimal curve and the initial guess (flip the graph upside-down on
%the y-axix); add a pretty legend and...profit!
hold on;
    plot(xOptPos,yOptPos,'ro-');
    set(gca,'YDir','reverse'); grid on;
    title('Brachistochrone','fontSize',14,'fontWeight','bold');
    xlabel('x','fontSize',14,'fontWeight','bold');
    ylabel('y   ','fontSize',14,'fontWeight','bold');
    set(get(gca,'YLabel'),'Rotation',0) %flip the axis
    %plotting the initial guess (curve).
    plot(xInitPos, yInitPos);
    legend( 'Optimal Curve', 'Initial Guess');
hold off;

%time step
dt = (tfOpt-t0)/N;

%Time domain used in this problem: N+1 points uniformly distributed on the
%interval [t0,tf] such that the first time node is t0 and the last sampled time
%node is tf
tDom = linspace(0,tfOpt,N+1);


%Plot the Lower Hamiltonian.  From Control theory, we know that the  Lower
%Hamiltonian is constant along the optimal curve, because this problem is time invarient,
%and from the hamiltonian value condition, we know that this constant is -1 (i.e.: we are
%minimizing final time).  Note that Fmul (of length numConstraints) is ordered as how F is
%packed, hence we start at index 2 since index 1 is the objective function tf.  Remenber,
%our dynamics were discritized as \dot x(t_i) = v(t_i)*sin(theta(T_i)) for i=0,...,N-1
%(same for yDot and vDot) Hence the straing looking indexing (cutting off the final time
%in the tDom)
Lx = Fmul(2:N+1);       %Legrange Mults for x-dynamics (N of them)
Ly = Fmul(N+2:2*N+1);   %Legrange Mults for y-dynamics (N of them)
Lv = Fmul(2*N+2:3*N+1); %Legrange Mults for v-dynamics (N of them)
DynX = vOpt(1:end-1).*sin(thetaOpt(1:end)); %time goes t_i for i=0,...,N-1
DynY = vOpt(1:end-1).*cos(thetaOpt(1:end)); %same as for x
DynV = (cos(thetaOpt(1:end)) - kFr*vOpt(1:end-1)); %same as for x
H = Lx.*DynX+Ly.*DynY+Lv.*DynV;
figure; plot(tDom(1:end-1),H,'*-')
title('Hamiltonian (should be approximatly -1)','fontSize',14,'fontWeight','bold');
xlabel('t','fontSize',14,'fontWeight','bold');
ylabel('H   ','fontSize',14,'fontWeight','bold');
set(get(gca,'YLabel'),'Rotation',0)
%--end plotting the Hamiltonian


%Plot the control (theta) which produces the optimal curve.  Quite valuable as
%sometimes the control can be non-smooth, and hence can cause problems when
%interpolating the control.  It's good practice to look at the control that
%produces the optimal solution.  Due to the discritization (Euler) process, theta was
%discritized to one less in size.
figure; plot( tDom(1:end-1), thetaOpt);
title('Theta (optimal) over time','fontSize',14,'fontWeight','bold');

%Save the variables for Feasibility Analysis on the Brach Problem
save( 'myVariables.mat',  ...
       'xOptPos', 'yOptPos', 'vOpt', 'thetaOpt', 'tfOpt', ...
       'tDom','N', 'x0','y0','v0','kFr');
   
%Call the script that performs feasibility analysis upon the solution returned
%by SNOPT only if Snopt found a solution that it liked (hence returning 1).  Note that
%sometimes Snopt fines solutions that it likes, but the solution is not the best solution.
if( INFO == 1 )   
    feasibilityAnalysis;
else
    disp('Feasibility Analysis skipped as INFO =/= 1');
end

%
% end brachFrictionJacobi.m
%