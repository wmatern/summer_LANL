function f = d_py(i,j,R)
order = 1;
if order == 1
% First order in the y direction
    f = (R(i,j+1)-R(i,j))./(R(i,j+1)+R(i,j)+eps);
elseif order == 2
% Second order in the y direction
    f = (-R(i,j+2)+4*R(i,j+1)-3*R(i,j))./(R(i,j+2)+4*R(i,j+1)+3*R(i,j)+eps)*2;
end
end