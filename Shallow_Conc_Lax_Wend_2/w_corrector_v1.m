function [ out ] = w_corrector(deltaT, dr, U_eigen, grad_half, grad_minus, grad_plus) 
%             deltaT       // Timestep
%             dr           // Cell's center to face distance
%             U_eigen      // State variable's eigenvalue (speed)
%             grad_half    // Centered gradient
%             grad_minus   // Downwind gradient
%             grad_plus    // Upwind gradient

   nu     = HALF * U_eigen * deltaT / dr;
   nu          = nu * (1.0 - nu);

   rdenom = ONE / max(SQR(grad_half), (real)EPSILON);
   rplus  = (grad_plus  * grad_half) * rdenom;
   rminus = (grad_minus * grad_half) * rdenom;

   out = .5*nu*(1.0- max(MIN3(1.0, rplus, rminus), 0.0));

end

