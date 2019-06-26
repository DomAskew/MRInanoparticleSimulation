function x = Diffusionwalk1 (n, dt, D, Bound, WalkX1 )
%
%  Stepsize is normal.
%
  s = sqrt(2.0 * D * dt)  * randn ( 1, n - 1 );
%
%  Compute the individual steps.
%
  x = WalkX1;
%
%  Direction is random.
%
dx(1,1:n-1) = s(1:n-1);
% 
%   Each position is the sum of the previous steps.
%            
% x(1,2:n) = cumsum ( dx(1,1:n-1) ); %original code before
% 
% Boundary condition.
%
    for dummy = 2:n
         x(1,dummy) = x(1,dummy-1) + dx(1,dummy-1);
       if (x(1,dummy)>Bound)
           x(1,dummy) = x(1,dummy)-(x(1,dummy)-Bound);
       elseif (x(1,dummy)<(-1*Bound))
           x(1,dummy) = x(1,dummy)-(x(1,dummy)+Bound);
       else
       end
    end
  return