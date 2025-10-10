function etaDot = eom(tau,eta,mu,dim3)
%% Purpose:
%
%  This code computes the equations of motion associated with the
%  circular restricted three-body problem (CR3BP) along with the
%  associated state transition matrix (STM) derivative with respect to
%  dimensionless time. This can be used with just N x 6 or N x 42 states
%  as long as the singleton dimension, dim3 is specified.
%
%% Inputs:
%
%  tau                      [N x 1]               Dimensionless Time
%
%  eta                      [N x 42]              Dimensionles State Vector
%                                                 [x,y,z,xdot,ydot,zdot,PHI]
%
%  mu                       double                Mass Ratio (m2/(m1+m2))
%
%  dim3                     integer               Singleton dimension
%                                                 specifier
%
%% Ouputs:
%
% etaDot                    [N x 42]              Derivative of the
%                                                 Dimensionless State
%                                                 Vector with respect to
%                                                 dimensionless time, tau
%
%
%% Revision History:
%  Darin C. Koblick                                              08/26/2025
%  Copyright 2025 Coorbital, Inc.
%% --------------------------- Begin Code Sequence ------------------------
if nargin == 0
              phi = eye(6);
               mu = 1.215058560962404E-2;
               x0 = [9.2294247990951594E-1	
                     0.0000000000000000E+0	
                     2.1592634764719765E-1	
                     4.2453615136710976E-13	
                     1.2786666755130519E-1	
                     -2.3343675990673871E-12
                     phi(:)];
             tau0 = 1.8036821222727220E+0;
             opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
          [tau,x] = ode113(@pumpkyn.cr3bp.eom,[0 tau0],x0,opts,mu,1);             
                M = reshape(x(end,7:42),[6 6]);    
               xf = x(end,1:6)';
          %compute the stability index:   
        lambdaMax = max(abs(eig(M)));
   stabilityIndex = 0.5.*(lambdaMax + 1./lambdaMax);
          figure('color',[1 1 1]);
          plot3(x(:,1),x(:,2),x(:,3));
          axis equal;
          title(['Stability Index = ',num2str(stabilityIndex)]);
          return;
end
%% Flatten the input dimension:
            [eta,fSeq] = pumpkyn.util.fDim(eta,dim3);
[x,y,z,xdot,ydot,zdot] = deal(eta(:,1),eta(:,2),eta(:,3), ...
                              eta(:,4),eta(:,5),eta(:,6));
%% Inputs: 
            N = size(eta,1);
           y2 = y.*y;
           z2 = z.*z;
           yz = y.*z;
            d = sqrt((x+mu).^2 + y2 + z2);
            r = sqrt((x-1+mu).^2 + y2 + z2);
           d3 = d.^3;
           r3 = r.^3;
       etaDot = eta.*0;

%% CR3B Equations of Motion:          
  etaDot(:,1) = xdot;
  etaDot(:,2) = ydot;
  etaDot(:,3) = zdot;
  etaDot(:,4) = -(1-mu).*(x+mu)./d3 - mu.*(x-1+mu)./r3 + 2.*ydot + x;
  etaDot(:,5) = -(1-mu).*y./d3 - mu.*y./r3 - 2.*xdot + y;
  etaDot(:,6) = -(1-mu).*z./d3 - mu.*z./r3;
  
%% Assemble the STM Derivative:
if size(eta,2) >= 42
    
     phi = reshape(permute(eta(:,7:42),[2 1]),[6 6 N]);
       F = phi.*0;
   term1 = (mu - 1)./((mu + x).^2 + y2 + z2).^1.5;
   term2 =  mu./((mu + x - 1).^2 + y2 + z2).^1.5;
   
F(4,1,:) = term1 - term2 + ...
           (3.*mu.*(2*mu + 2.*x - 2).*(mu + x - 1))./ ...
           (2.*((mu + x - 1).^2 + y2 + z2).^2.5) - ...
           (3.*(2.*mu + 2.*x).*(mu + x).*(mu - 1))./ ...
           (2.*((mu + x).^2 + y2 + z2).^2.5) + 1; 
     
F(5,1,:) = (3.*mu.*y.*(mu + x - 1))./((mu + x - 1).^2 + y2 + z2).^2.5 - ...
           (3.*y.*(mu + x).*(mu - 1))./((mu + x).^2 + y2 + z2).^2.5;     
       
F(6,1,:) = (3.*mu.*z.*(mu + x - 1))./((mu + x - 1).^2 + y2 + z2).^2.5 - ...
           (3.*z.*(mu + x).*(mu - 1))./((mu + x).^2 + y2 + z2).^2.5;  
     
F(4,2,:) = F(5,1,:);
 
F(5,2,:) = term1 - term2 - ...
          (3.*y2.*(mu - 1))./((mu + x).^2 + y2 + z2).^2.5 + ...
          (3.*mu.*y2)./((mu + x - 1).^2 + y2 + z2).^2.5 + 1;   
           
F(6,2,:) = (3.*mu.*yz)./((mu + x - 1).^2 + y2 + z2).^2.5 - ...
           (3.*yz*(mu - 1))./((mu + x).^2 + y2 + z2).^2.5; 
     
F(4,3,:) = F(6,1,:);
     
F(5,3,:) = F(6,2,:); 
      
F(6,3,:) = term1 - term2 - ...
          (3.*z2.*(mu - 1))./((mu + x).^2 + y2 + z2).^2.5 + ...
          (3.*mu.*z2)./((mu + x - 1).^2 + y2 + z2).^2.5;
      
  F(1:3,4:6,:) = repmat(eye(3),[1 1 N]);
  F(4:6,4:6,:) = repmat([0 2 0; -2 0 0; 0 0 0],[1 1 N]);
  if size(F,3) == 1 && size(phi,3) == 1
        phiDot = F*phi;   
  else
        phiDot = mntimes(F,phi,1,2,1,2); %F*phi;   
  end
etaDot(:,7:42) = permute(pumpkyn.util.fDim(phiDot,3),[2 1]); 
end

%% Reshape to original dimensions:
etaDot = pumpkyn.util.eDim(etaDot,fSeq);
end