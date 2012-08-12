function f = d_px(j,i,R)
order = 1;
if order == 1
% First Order in the x direction
    f = (R(j,i+1)-R(j,i  ))./(R(j,i+1)+R(j,i  )+eps);
elseif order == 2
% Second order in the x direction
    f = (-R(j,i+2)+4*R(j,i+1)-3*R(j,i))./(R(j,i+2)+4*R(j,i+1)+3*R(j,i)+eps)*2;
end
end