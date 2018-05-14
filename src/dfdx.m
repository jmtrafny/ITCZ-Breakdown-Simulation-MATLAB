function dfdx = x_deriv_func(phi, dx)
%Takes derivative with respect to x

[m,n] = size(phi); %Check dimensions of data
dfdx(:,1) = (phi(:,2)-phi(:,n))/(2*dx); %Backward diff for last value
dfdx(:,n) = (phi(:,1)-phi(:,n-1))/(2*dx); %Forward diff for first value
dfdx(:,2:n-1) = (phi(:,3:n)-phi(:,1:n-2))/(2*dx); %Center diff

end

