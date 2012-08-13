
function ddx = d_dx(phi,i,j)
    ddx = 1/12*phi(i,j-2)-2/3*phi(i,j-1)+2/3*phi(i,j+1)-1/12*phi(i,j+2);
end