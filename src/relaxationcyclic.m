function [a] = relaxationcyclic(GN,zeta,dx,tol)
%{
Function that calculates the solution to Poisson's equation
Laplacian(a)=zeta where a is the streamfunction. It outputs the stream
function (psi) at the new time step given the vorticity at the new time
step and a first guess for the stream function (i.e., psi at the current
time step). Inputs are the psi field, zeta, grid space size, and tol (the
accuracy tolerance or maximum error allowed. 

Credit to Tyler Green for developing the inital version of this code. 
%}

[imax,jmax] = size(GN);
[imax2,jmax2] = size(zeta);

if ((imax ~= imax2) || (jmax ~=jmax2))
    warning('Array dimensions do not match. Exiting')
    return
end

t = cos(pi/jmax)+cos(pi/imax);
alpha = (8-4*sqrt(4-t^2))/t^2;
%  alpha =1;



%% Iterate on initial guess until tolerance is met

residual = zeros(imax,jmax);
residual(3:imax-2,1:jmax) = 1;

 while max(max(abs(residual))) >tol

    for i = 2:imax-1
        for j = 1:jmax
            %Calculate residual and update guess at first column using cyclic
            %boundary conditions 
            if j==1
                residual(i,j)= ((GN(i+1,j)+GN(i-1,j)+GN(i,j+1)+GN(i,jmax)-4*GN(i,j))/dx^2)-zeta(i,j);
                GN(i,j) = GN(i,j)+alpha*(dx^2/4)*residual(i,j);
            %Calculate residual and update guess at interior points
            elseif j>1 && j<jmax
                residual(i,j) = ((GN(i+1,j)+GN(i-1,j)+GN(i,j+1)+GN(i,j-1)-4*GN(i,j))/dx^2)-zeta(i,j);             
                GN(i,j) = GN(i,j)+alpha*(dx^2/4)*residual(i,j);
            %Calculate residual and update guess at last column using cyclic
            %boundary conditions 
            else 
                residual(i,j) = ((GN(i+1,j)+GN(i-1,j)+GN(i,1)+GN(i,jmax-1)-4*GN(i,j))/dx^2)-zeta(i,j);
                GN(i,j) = GN(i,j)+alpha*(dx^2/4)*residual(i,j);
            end
        end
    end
   
    
end

%% Assign GN to a, the outputted streamfunction
a = GN;
% fprintf('It took %d iterations to converge\n',n)
end