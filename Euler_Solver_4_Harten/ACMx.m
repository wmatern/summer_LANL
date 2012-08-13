function w = ACMx(u,i,j,a)

alpha(j,i) = a*(abs(u(j,i+1)-2*u(j,i)+u(j,i-1))./...
    (abs(u(j,i+1)-u(j,i))+abs(u(j,i)-u(j,i-1))+eps)).^2;

w = minmod(alpha(j,i)/2.*(minmod(u(j,i+1)+u(j,i),u(j,i-1)+u(j,i))),...
    minmod(u(j,i+1)-(u(j,i+1)+u(j,i))/2,u(j-1,i)-(u(j,i-1)+u(j,i))/2));