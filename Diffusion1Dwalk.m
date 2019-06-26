function x = Diffusion1Dwalk ( m, dt, d, t, bound )
%
%  Number of steps, n:
%
  n =( t / dt)+1;
%
%  Compute the individual steps.
%
  x = zeros ( m, n );
%
%  Stepsize is normal.
%
  s =  sqrt(2.0 * d * dt)  * randn ( 1, n - 1 );
%
%  Direction is random.
%
dx(1:m,1:n-1) = s(1:n-1);
%   if ( m == 1 )
%     dx(1:m,1:n-1) = s(1:n-1);
%   else
%     a = randn ( m, n - 1 );
%     v = s ./ sqrt ( sum ( a.^2 ) );
%     b = spdiags ( v', 0, n-1, n-1 );
%     dx(1:m,1:n-1) = a * b;
%   end
%
%  Each position is the sum of the previous steps.
%
 dummy = 2;
    for dummy = 2:n
         x(1:m,dummy) = x(1:m,dummy-1) + dx(1:m,dummy-1);
       if (x(1:m,dummy)>bound)
           x(1:m,dummy) = x(1:m,dummy)-(x(1:m,dummy)-bound);
       elseif (x(1:m,dummy)<(-1*bound))
           x(1:m,dummy) = x(1:m,dummy)-(x(1:m,dummy)+bound);
       else
       end
    end
%                
%   x(1:m,2:n) = cumsum ( dx(1:m,1:n-1), 2 ); %original code before
%    boundary conditions.

  return