
function ddy = d_dy(phi,i,j)
        ddy = 1/12*phi(i-2,j)-2/3*phi(i-1,j)+2/3*phi(i+1,j)-1/12*phi(i+2,j);
end