%% Purpose:
%
%  This routine will symbolically compute the jacobian matrix ( DyF(y(t)) ) 
%  corresonding to the equations of motion present in tfMin.
%

% Clear everyhthing:
clear all;

%Declare initial variables:
syms x y z xDot yDot zDot lambda_rx lambda_ry lambda_rz real;
syms lambda_vx lambda_vy lambda_vz lambda_m muStar m c u Tmax real;

%Define costates:
lambda_r = [lambda_rx; lambda_ry; lambda_rz];
lambda_v = [lambda_vx; lambda_vy; lambda_vz];
    rDot = [xDot; yDot; zDot];

%Define g(r) and h(v) functions:
       d = sqrt((x + muStar)^2 + y^2 + z^2);
       r = sqrt((x - 1 + muStar)^2 + y^2 + z^2);
      d3 = d^3;
      r3 = r^3;
      d5 = d^5;
      r5 = r^5;
      r7 = r^7;
      d7 = d^7;
      gr = [x - (1-muStar)*(x+muStar)/d3 - muStar*(x-1+muStar)/r3
            y - (1-muStar)*y/d3 - muStar*y/r3
            0 - (1-muStar)*z/d3 - muStar*z/r3];
      hv = [+2.*yDot
            -2.*xDot;
            +0];

%Determine G = d g(r) / dr
      G = jacobian(gr,[x,y,z]);

%Determine H = d h(v) / dv
      H = jacobian(hv,[xDot,yDot,zDot]); 

%Define the state derivative function F(y)
       alpha = -lambda_v./norm(lambda_v);                  %Equation 10
        vDot =  gr + hv + u.*alpha.*Tmax./m;
        mDot = -u*Tmax/c;

%Costate Derivatives:
 lambdaDot_r = -G'*lambda_v;
 lambdaDot_v = -lambda_r - H'*lambda_v;
 lambdaDot_m = -norm(lambda_v)*u*Tmax/m^2;
 
 %All State Derivatives (Equation 15):
           Fy = [rDot(:); 
                 vDot(:); 
                 mDot; 
                 lambdaDot_r(:); 
                 lambdaDot_v(:); 
                 lambdaDot_m];

 %Compute the Derivative of the State Transition Matrix (STM):
        DyF =  jacobian(Fy, [x,y,z,xDot,yDot,zDot,m,lambda_r(:)',lambda_v(:)',lambda_m]);   
        DyF =  simplify(DyF,'IgnoreAnalyticConstraints', true);
        
 %Substittue easy expressions, then simplify:       
        clear r3 r5 r7 d3 d5 d7;
        DyF = subs(DyF,((muStar + x - 1)^2 + y^2 + z^2)^(3/2),'r3');
        DyF = subs(DyF,((muStar + x - 1)^2 + y^2 + z^2)^(5/2),'r5');
        DyF = subs(DyF,((muStar + x - 1)^2 + y^2 + z^2)^(7/2),'r7');
        DyF = subs(DyF,((muStar + x)^2 + y^2 + z^2)^(3/2),'d3');
        DyF = subs(DyF,((muStar + x)^2 + y^2 + z^2)^(5/2),'d5');
        DyF = subs(DyF,((muStar + x)^2 + y^2 + z^2)^(7/2),'d7');
        DyF = subs(DyF,(lambda_vx^2 + lambda_vy^2 + lambda_vz^2)^(1/2),'lambda_v_Mag');
        DyF =  expand(DyF);
        DyF =  simplify(DyF,'IgnoreAnalyticConstraints', true);

