%{
    The sparsity pattern of G: The Jacobian of the constraints (packed into F)
    with respect to the Decesion Variables (packed into xInit) tells Snopt which
    sections of the jacobian that are nonZero, and hence need to be computed
    (either by me, or by Snopt; in this version of the code I explicitly compute
    the entire Jacobian).  That is, the Jacobian is usually a giant matrix; just
    for this small Brach problem, the Full Jacobian is 

    Note: As mentioned in the paper, (and in the main script, brachFriction.m)
          the order of the Jacobian and it's sparsity is determined by the order
          of how
            * The decision variables were packed into xInit (in
              buildInitialGuess)
            * The constraint variables were packed into F (in userFun)
    Inputs:  
        D: The Euler Diff Matrix, which is the sparsity pattern of sections of
           the Jacobian: D(dynamicsX) wrt x, D(dynamicsY) wrt y, and 
           D(dynamicsV) wrt v
        N: the number of sample points (N+1) on [t0,tf]; just knowing N we can
           build the sparsity pattern of the Jacobian.
    Outputs: A nD by nC matrix M (where nD is the number of decesion variables
             and nC is the number of constraints (packed into F) such that
                M(i,j) =/= 0 iff this is an element of the Jacobian that is not
                                 zero
             Hence, applying find on M returns the locations of the nonzero
             elements in M.
%}
function [ jacobianSparsity ] = findSparsityPattern( D, N )

    %S is the sparsity for matrices alpha and delta
    %A and B are the sparsity patterns for the initial and final conditions on x
    %and y  eg 1 by N+1 row vector A==[10...0] B==[0...01]
    S = horzcat( eye(N), zeros(N,1));
    A = horzcat( 1, zeros(1,N));
    B = horzcat( zeros(1,N), 1);
    
    
    %Directly building the sparsity patern (a 1 where we need to calculate the
    %jacobian) is numConstraints by numDecVars. This sparsity follows the
    %structure given in the Brach paper.
    jacobianSparsity = ...
        [ zeros(1, N+1) zeros(1, N+1) zeros(1, N+1) zeros(1, N)    1;      ...
              D         zeros(N, N+1)     S           eye(N)      ones(N,1); ...
          zeros(N, N+1)     D             S           eye(N)      ones(N,1); ...
          zeros(N, N+1) zeros(N, N+1)     D           eye(N)      ones(N,1); ...
              A         zeros(1, N+1) zeros(1, N+1) zeros(1, N)    0;      ...
          zeros(1, N+1)     A         zeros(1, N+1) zeros(1, N)    0;      ...
          zeros(1, N+1) zeros(1, N+1)     A         zeros(1, N)    0;      ...
              B         zeros(1, N+1) zeros(1, N+1) zeros(1, N)    0;      ...
          zeros(1, N+1)     B         zeros(1, N+1) zeros(1, N)    0 ];

end

