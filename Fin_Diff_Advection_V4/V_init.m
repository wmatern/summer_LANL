function [ V ] = V_init( X,Y,nx,ny )
V = zeros(nx+4,ny+4);
xnorm = reshape(X,nx,ny)/nx;
ynorm = reshape(Y,nx,ny)/ny;

V(3:nx+2,3:ny+2) = -(10*sqrt((xnorm-.5).^2 + (ynorm-.5).^2)) .* (-(xnorm-.5)./sqrt((xnorm-.5).^2+(ynorm-.5).^2));
%V(3:nx+2,3:ny+2) = 0;

V (:,2  ) = V (:,3  );
V (:,ny+3) = V (:,ny+2);
V (2  ,:) = V (3  ,:);
V (nx+3,:) = V (nx+2,:);

V (:,1  ) = V (:,2  );
V (:,ny+4) = V (:,ny+3);
V (1  ,:) = V (2  ,:);
V (nx+4,:) = V (nx+3,:);
end

