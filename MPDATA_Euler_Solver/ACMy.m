function w = ACMy(u,i,j,a)

alpha(j,i) = a*(abs(u(j+1,i)-2*u(j,i)+u(j-1,i))./...
    (abs(u(j+1,i)-u(j,i))+abs(u(j,i)-u(j-1,i))+eps)).^2;

w = minmod(alpha(j,i)/2.*(minmod(u(j+1,i)+u(j,i),u(j-1,i)+u(j,i))),...
    minmod(u(j+1,i)-(u(j+1,i)+u(j,i))/2,u(j-1,i)-(u(j-1,i)+u(j,i))/2));

% u_L(j,i) = u(j,i) - 