function [ U ] = U_init( X,Y,nx,ny )
U = zeros(nx+4,ny+4);
xnorm = reshape(X,nx,ny)/nx;
ynorm = reshape(Y,nx,ny)/ny;

U(3:nx+2,3:ny+2) = (10*sqrt((xnorm-.5).^2 + (ynorm-.5).^2)) .* ((ynorm-.5)./sqrt((xnorm-.5).^2+(ynorm-.5).^2));
%U(3:nx+2,3:ny+2) = 1;

U (:,2  ) = U (:,3  );
U (:,ny+3) = U (:,ny+2);
U (2  ,:) = U (3  ,:);
U (nx+3,:) = U (nx+2,:);

U (:,1  ) = U (:,2  );
U (:,ny+4) = U (:,ny+3);
U (1  ,:) = U (2  ,:);
U (nx+4,:) = U (nx+3,:);
end

