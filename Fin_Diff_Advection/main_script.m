clear all; clc; close all

% Model Parameters
a     = 2;                % slope modification constant 1.9<a<2.3
kappa = 1/3;              % MUSCL interpolant parameter

% Switches
Slope_Mod = 1;            % 1 means on
BC        = 2;            % 1=reflection, 2=periodic
ic        = 2;            % 1=zero vel, 2=only u vel, 3=rotating
SCH       = 1;            % 1=LaxWnd, 2=MPDATA, 3=MUSCL

% Problem Parameters
nx = 64;                  % grid size
ny = 64;                  % grid size
n  = 64;                  % FIX THIS LATER!!!!!!!!!!!!!!
g  = 9.8;                 % gravitational constant
dt = 0.02;                % hardwired timestep
dx = 1.0;
dy = 1.0;
Cx = dt/dx;
Cy = dt/dy;

nplotstep = 8;            % plot interval

% Initialize graphics
% [surfplot,phiplot,heightplot,top] = initgraphics (n,nx,ny);

% Outer loop, restarts.

if ic == 1
    H = ones (n+4,n+4); U = zeros (n+4,n+4); V = zeros (n+4,n+4);
elseif ic == 2
    H = ones (n+4,n+4); U = ones  (n+4,n+4); V = zeros (n+4,n+4);
elseif ic == 3
    H = ones (n+4,n+4); U = ones (n+4,n+4); V = zeros (n+4,n+4);
    U (3:n+2,3:n+2) = sqrt(X.^2./Y.^2); V (3:n+2,3:n+2) = Y./X;
end
phi          = zeros (n+4,n+4);
phi(1:n/2+2,:) = 1;
nstep = 0;

i = 3:n+2; j = 3:n+2;
plot1(1:length(j)  ,1) = phi(j,34);
if BC == 2, plot1(length(j)+j-2,1) = phi(j,34); end
figure(1),plot(plot1,'.-') , title('phi')

% Main loop
while nstep < 10000000
    nstep = nstep + 1;
    
    if BC == 1 % Boundary condition switch
        phi (:,2  ) = phi (:,3  );
        phi (:,n+3) = phi (:,n+2);
        phi (2  ,:) = phi (3  ,:);
        phi (n+3,:) = phi (n+2,:);
        
        phi (:,1  ) = phi (:,2  );
        phi (:,n+4) = phi (:,n+3);
        phi (1  ,:) = phi (2  ,:);
        phi (n+4,:) = phi (n+3,:);
    elseif BC == 2
        % Periodic boundary conditions
        phi = Periodic_BCx(phi, nx);
        phi = Periodic_BCy(phi, ny);
    end
    
    if SCH == 1
        % First half step
       % x direction
       i = 2:n+2;
       j = 2:n+1;
       
       % phi
       phix(i,j) = (phi(i+1,j+1) + phi(i,j+1))/2 - ...
           dt/(4*dx) * (U(i+1,j+1)./H(i+1,j+1)+U(i,j+1)./H(i,j+1)) .*...
           (phi(i+1,j+1)-phi(i,j+1));
       
       % height
       Hx (i,j) = (H (i+1,j+1)+H (i,j+1))/2 ...
           - dt/(2*dx)*(U (i+1,j+1)-U (i,j+1));
   
       % x momentum
       Ux (i,j) = (U (i+1,j+1)+U (i,j+1))/2 -  ...
           dt/(2*dx)*((U (i+1,j+1).^2./H(i+1,j+1) + g/2*H (i+1,j+1).^2) - ...
           (U (i,j+1).^2./H(i,j+1) + g/2*H (i,j+1).^2));
       
       % y momentum
       Vx (i,j) = (V (i+1,j+1)+V (i,j+1))/2 - ...
           dt/(2*dx)*((U (i+1,j+1).*V (i+1,j+1)./H(i+1,j+1)) - ...
           (U (i,j+1).*V (i,j+1)./H(i,j+1)));
       
       % y direction
       i = 2:n+1;
       j = 2:n+2;
      
       % phi
       phiy(i,j) = (phi(i+1,j+1) + phi(i+1,j))/2 -...
           dt/(4*dy)*(V(i+1,j+1)./H(i+1,j+1)+V(i+1,j)./H(i+1,j)) .*...
           (phi(i+1,j+1)-phi(i+1,j));
       
       % height
       Hy (i,j) = (H (i+1,j+1)+H (i+1,j))/2 - dt/(2*dy)*(V (i+1,j+1)-V (i+1,j));
       
       % x momentum
       Uy (i,j) = (U (i+1,j+1)+U (i+1,j))/2 - ...
           dt/(2*dy)*((V (i+1,j+1).*U (i+1,j+1)./H(i+1,j+1)) - ...
           (V (i+1,j).*U (i+1,j)./H(i+1,j)));
       % y momentum
       Vy (i,j) = (V (i+1,j+1)+V (i+1,j))/2 - ...
           dt/(2*dy)*((V (i+1,j+1).^2./H(i+1,j+1) + g/2*H (i+1,j+1).^2) - ...
           (V (i+1,j).^2./H(i+1,j) + g/2*H (i+1,j).^2));
   
       % Second half step
       i = 3:n+2;
       j = 3:n+2;
       
       
       % concentration
       nuplusy  = abs(Vy(i-1,j  )./Hy(i-1,j  ))*Cy;
       nuplusx  = abs(Ux(i  ,j-1)./Hx(i  ,j-1))*Cx;
       numinusy = abs(Vy(i-1,j-1)./Hy(i-1,j-1))*Cy;
       numinusx = abs(Ux(i-1,j-1)./Hx(i-1,j-1))*Cx;

       phi (i,j) = phi (i,j) - (dt/(2*dx))*(Ux (i,j-1)./Hx(i,j-1)...
           + Ux (i-1,j-1)./Hx(i-1,j-1)) .* (phix(i,j-1)-phix(i-1,j-1))...
           - (dt/(2*dy))*(Vy (i-1,j)./Hy(i-1,j)...
           + Vy (i-1,j-1)./Hy(i-1,j-1)) .* (phiy(i-1,j)-phiy(i-1,j-1))...
           - wminusx(phi(:,:),i,j,numinusx,2).*(phi(i  ,j  )-phi(i-1,j  ))...
           + wplusx (phi(:,:),i,j,nuplusx ,2).*(phi(i+1,j  )-phi(i  ,j  ))...
           - wminusy(phi(:,:),i,j,numinusy,2).*(phi(i  ,j  )-phi(i  ,j-1))...
           + wplusy (phi(:,:),i,j,nuplusy ,2).*(phi(i  ,j+1)-phi(i  ,j  ));
    elseif SCH == 2
        U_halfp (i,j) = (U(i,j+1)./H(i,j+1) + U(i,j  )./H(i,j  ))*Cx/2;
        U_halfm (i,j) = (U(i,j  )./H(i,j  ) + U(i,j-1)./H(i,j-1))*Cx/2;
        V_halfp (i,j) = (V(i+1,j)./H(i+1,j) + V(i,j  )./H(i,j  ))*Cy/2;
        V_halfm (i,j) = (V(i,j  )./H(i,j  ) + V(i-1,j)./H(i-1,j))*Cy/2;
        
        phi_New (i,j) = phi(i,j) - (MPDATA(phi(i,j  ),phi(i,j+1),V_halfp(i,j))...
                                  - MPDATA(phi(i,j-1),phi(i,j  ),V_halfm(i,j)))...
                                 - (MPDATA(phi(i,j  ),phi(i+1,j),U_halfp(i,j))...
                                  - MPDATA(phi(i-1,j),phi(i,j  ),U_halfm(i,j)));
        
        if BC == 2
            phi_New = Periodic_BCx(phi_New, nx);
            phi_New = Periodic_BCy(phi_New, ny);
        end
        
        U_d_p   (i,j) = dt/dx*((abs(U_halfp(i,j))-U_halfp(i,j).^2).*...
            d_px(i,j  ,phi_New) - U_halfp(i,j).*V_halfp(i,j).*d_py(i,j  ,phi_New));
        U_d_m   (i,j) = dt/dx*((abs(U_halfm(i,j))-U_halfm(i,j).^2).*...
            d_px(i,j-1,phi_New) - U_halfm(i,j).*V_halfm(i,j).*d_py(i,j-1,phi_New));
        V_d_p   (i,j) = dt/dy*((abs(V_halfp(i,j))-V_halfp(i,j).^2).*...
            d_py(i  ,j,phi_New) - U_halfp(i,j).*V_halfp(i,j).*d_px(i  ,j,phi_New));
        V_d_m   (i,j) = dt/dy*((abs(V_halfm(i,j))-V_halfm(i,j).^2).*...
            d_py(i-1,j,phi_New) - U_halfm(i,j).*V_halfm(i,j).*d_px(i-1,j,phi_New));
        
        phi (i,j) = phi_New(i,j)...
            - (MPDATA(phi_New(i,j  ),phi_New(i,j+1),V_d_p(i,j))...
            -  MPDATA(phi_New(i,j-1),phi_New(i,j  ),V_d_m(i,j)))...
            - (MPDATA(phi_New(i,j  ),phi_New(i+1,j),U_d_p(i,j))...
            -  MPDATA(phi_New(i-1,j),phi_New(i,j  ),U_d_m(i,j)));
    elseif SCH == 3
        % phi
        Uxm = (U(i  ,j  )./H(i  ,j  ) + U(i-1,j  )./H(i-1,j  ))/2;
        Vym = (V(i  ,j  )./H(i  ,j  ) + V(i  ,j-1)./H(i  ,j-1))/2;
        Uxp = (U(i+1,j  )./H(i+1,j  ) + U(i  ,j  )./H(i  ,j  ))/2;
        Vyp = (V(i  ,j+1)./H(i  ,j+1) + V(i  ,j  )./H(i  ,j  ))/2;
        
        
        % X direction
        phixpL1 = phi(i  ,j) + 0.25*((1+kappa)*(phi(i+1,j)...
            - phi(i  ,j))+(1-kappa)*(phi(i  ,j)-phi(i-1,j)));
        phixpR1 = phi(i+1,j) - 0.25*((1+kappa)*(phi(i+1,j)...
            - phi(i  ,j))+(1-kappa)*(phi(i+2,j)-phi(i+1,j)));
        phixmL1 = phi(i-1,j) + 0.25*((1+kappa)*(phi(i  ,j)...
            - phi(i-1,j))+(1-kappa)*(phi(i-1,j)-phi(i-2,j)));
        phixmR1 = phi(i  ,j) - 0.25*((1+kappa)*(phi(i  ,j)...
            - phi(i-1,j))+(1-kappa)*(phi(i+1,j)-phi(i  ,j)));
        
        phixpL = phi(i,j) + minmod3(phixpL1 - phi(i,j),...
            phi(i+1,j)-phi(i,j),phi(i,j)-phi(i-1,j));
        phixpR = phi(i+1,j) - minmod3(phi(i+1,j)-phixpR1,...
            phi(i+1,j)-phi(i,j), phi(i+2,j)-phi(i+1,j));
        phixmL = phi(i-1,j) + minmod3(phixmL1 - phi(i-1,j),...
            phi(i,j)-phi(i-1,j),phi(i-1,j)-phi(i-2,j));
        phixmR = phi(i,j) - minmod3(phi(i,j)-phixmR1,...
            phi(i,j)-phi(i-1,j), phi(i+1,j)-phi(i,j));
        
        if Slope_Mod == 2 % Slope modification switch
            dphix = minmod(a*minmod(phixmR-phixmL,phixpR-phixpL),...
                minmod(phi(i+1,j)-phixpL,phixmR-phi(i-1,j)));
        else
            dphix = 0;
        end
        
        phixpL = phixpL + dphix;
        phixpR = phixpR - dphix;
        phixmL = phixmL + dphix;
        phixmR = phixmR - dphix;
        
        fxpL = Uxp.*phixpL;
        fxpR = Uxm.*phixpR;
        fxmL = Uxp.*phixmL;
        fxmR = Uxm.*phixmR;
        
        betap = max(max(Uxp, Uxp), eps);
        betam = max(max(Uxm, Uxm), eps);
        
        fxp = 0.5*(fxpL + fxpR  - betap.*(phixpR-phixpL));
        fxm = 0.5*(fxmL + fxmR  - betam.*(phixmR-phixmL));
        
        % Y direction
        phiypL1 = phi(i,j  ) + 0.25*((1+kappa)*(phi(i,j+1)...
            - phi(i,j  ))+(1-kappa)*(phi(i,j  )-phi(i,j-1)));
        phiypR1 = phi(i,j+1) - 0.25*((1+kappa)*(phi(i,j+1)...
            - phi(i,j  ))+(1-kappa)*(phi(i,j+2)-phi(i,j+1)));
        phiymL1 = phi(i,j-1) + 0.25*((1+kappa)*(phi(i,j  )...
            - phi(i,j-1))+(1-kappa)*(phi(i,j-1)-phi(i,j-2)));
        phiymR1 = phi(i,j  ) - 0.25*((1+kappa)*(phi(i,j  )...
            - phi(i,j-1))+(1-kappa)*(phi(i,j+1)-phi(i,j  )));
        
        phiypL = phi(i,j) + minmod3(phiypL1 - phi(i,j),...
            phi(i,j+1)-phi(i,j),phi(i,j)-phi(i,j-1));
        phiypR = phi(i,j+1) - minmod3(phi(i,j+1)-phiypR1,...
            phi(i,j+1)-phi(i,j), phi(i,j+2)-phi(i,j+1));
        phiymL = phi(i,j-1) + minmod3(phiymL1 - phi(i,j-1),...
            phi(i,j)-phi(i,j-1),phi(i,j-1)-phi(i,j-2));
        phiymR = phi(i,j) - minmod3(phi(i,j)-phiymR1,...
            phi(i,j)-phi(i,j-1), phi(i,j+1)-phi(i,j));
        
        if Slope_Mod % Slope Modification switch
            dphiy = minmod(a*minmod(phiymR-phiymL,phiypR-phiypL),...
                minmod(phi(i,j+1)-phiypL,phiymR-phi(i,j-1)));
        else
            dphiy = 0;
        end
        
        phiypL = phiypL + dphiy;
        phiypR = phiypR - dphiy;
        phiymL = phiymL + dphiy;
        phiymR = phiymR - dphiy;
        
        fypL = Vyp.*phiypL;
        fypR = Vym.*phiypR;
        fymL = Vyp.*phiymL;
        fymR = Vym.*phiymR;
        
        betap = max(max(Vyp, Vyp), eps);
        betam = max(max(Vym, Vym), eps);
        
        fyp = 0.5*(fypL + fypR  - betap.*(phiypR-phiypL));
        fym = 0.5*(fymL + fymR  - betam.*(phiymR-phiymL));
        
        % Integrate concentration
        phi(i,j) = phi(i,j) + Cx*(fxm - fxp)...
            + Cy*(fym - fyp);
    end
        
    % Update plot
    if mod (nstep,nplotstep) == 0
        plot1(1:length(j)  ,1) = phi(j,34);
        if BC == 2, plot1(length(j)+j-2,1) = phi(j,34); end

        figure(1),plot(plot1,'.-') , title('phi')
    end
%         if any (any (isnan (H))), break, end  % Unstable, restart
%         if any (any (isinf (H))), break, end  % Unstable, restart
%         if nstep>80, break, end
end