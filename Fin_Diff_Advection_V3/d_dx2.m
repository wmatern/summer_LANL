
function ddx = d_dx2(phi,i,j)
    r = (phi(i,j)-phi(i,j-1))/(phi(i,j+1)-phi(i,j));
    fee = max(0,min(1,r));
    
    f_low_minus = f_low(phi,i,j-1);
    f_low_plus  = f_low(phi,i,j  );
    
    f_high_minus = f_high(phi,i,j-1);
    f_high_plus  = f_high(phi,i,j  );
    
    F_plus  = f_low-fee*(f_low+f_high);
    F_minus = f_low
end

function f = f_low(phi,i,j)
    f = (phi(i,j+1)+phi(i,j))/2;
end

function f = f_high(phi,i,j)
    f = 1/12*phi(i,j-2)-2/3*phi(i,j-1)+2/3*phi(i,j+1)-1/12*phi(i,j+2);
end