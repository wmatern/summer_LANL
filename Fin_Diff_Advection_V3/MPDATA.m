function F = MPDATA(R_L, R_R, U)
Up = 0.5*(U+abs(U));
Um = 0.5*(U-abs(U));

F = Up.*R_L+Um.*R_R;
end