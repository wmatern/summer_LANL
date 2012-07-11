function f = d_py(j,i,R)
order = 2;
if order == 1
% First order in the y direction
    f = (R(j+1,i)-R(j  ,i))./(R(j+1,i)+R(j  ,i)+eps);
elseif order == 2
% Second order in the y direction
    f = (-R(j+2,i)+4*R(j+1,i)-3*R(j,i))./(R(j+2,i)+4*R(j+1,i)+3*R(j,i)+eps)*2;
end
end