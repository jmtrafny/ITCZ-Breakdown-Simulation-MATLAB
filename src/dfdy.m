function dfdy = y_deriv_func(phi, dy)
%Takes derivative with respect to x

[m,n] = size(phi); %Check dimensions of data
dfdy=zeros(m,n);
dfdy(2:m-1,:) = (phi(3:m,:)-phi(1:m-2,:))/(2*dy); %Center diff

end

