function RHO = Periodic_BCy(RHO, ny)
RHO(2,:) = RHO(ny+2,:);
RHO(1,:) = RHO(ny+1,:);
RHO(ny+3,:) = RHO(3,:);
RHO(ny+4,:) = RHO(4,:);
end