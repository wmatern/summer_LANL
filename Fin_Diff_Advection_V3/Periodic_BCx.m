function RHO = Periodic_BCx(RHO, nx)
RHO(:,2) = RHO(:,nx+2);
RHO(:,1) = RHO(:,nx+1);
RHO(:,nx+3) = RHO(:,3);
RHO(:,nx+4) = RHO(:,4);
end