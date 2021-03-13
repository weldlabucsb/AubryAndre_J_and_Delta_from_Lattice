function [J,Delta] = J_Delta_PiecewiseFit(s1,s2)
%A piecewise fit of the numerically computed values of J and Delta (from
%J_Numeric_Approx.m) - only applicapable for lambda1=1064 and lambda2=915.
%
%   s1 is the 1064nm laser lattice depth in units of the recoil energy of a
%   photon with wavelength 1064nm.  Note that this must be the primary
%   lattice
%   s2 is the 915nm laser lattice depth in units of the recoil energy of the 
%   recoil energy of a photon with wavelength 915nm
%
%   The max error among all of the fits is 0.1% compared to the numerically
%   calculated values.  The fit functions only depend on s1 since Delta is
%   otherwise just proportional to s2.
%
%   See Michele Modugno 2009 New J. Phys. 11, 033023 for details of the
%   theory.  For Delta, only the integral over wannier states is determined
%   by the fit.

assert((length(s1)==1)&(length(s2)==1), 'This function does not support vector inputs');

la1 = 1064;
la2 = 915;

if ((s1>=17)&&(s1<=40))
    % s1 between 17Er and 40Er

    % Fit parameters for J
    A = 1.907777300434345;
    alpha = 1.110634053285819;
    B = 2.378458268375515;
    beta = 0.478324506582921;
    
    % Calculating J
    J = A*s1^alpha*exp(-B*s1^beta); % Error less than 0.011% in this range.
    
    
    % Fit parameters for Delta
    C = 0.967975012322732;
    D = 2.062460149269072;
    gamma = 0.631843655943315;
    
    % Calculating Delta
    Delta = 0.5*(s2*(la1/la2)^2)*C*exp(-D*s1^(-gamma)); % Error less than 0.011% in this range.
    
elseif ((s1>=12)&&(s1<17))
    % s1 between 12Er and 17Er
    
    % Fit parameters for J
    A = 1.773300418805690;
    alpha = 1.128356804009188;
    B = 2.351229520668319;
    beta = 0.481494822456852;
    
    % Calculating J
    J = A*s1^alpha*exp(-B*s1^beta);  % Error less than 0.05% in this range.
    
    
    % Fit parameters for Delta
    C = 0.928133240959474;
    D = 2.326438686265068;
    gamma = 0.720086338581454;
    
    % Calculating Delta
    Delta = 0.5*(s2*(la1/la2)^2)*C*exp(-D*s1^(-gamma)); % Error less than 0.011% in this range.
    
elseif ((s1>=8)&&(s1<12))
    % s1 between 8Er and 12Er
    
        
    % Fit parameters for J
    A = 3.315925688116741;
    alpha = 1.469663378695678;
    B = 3.088737434881095;
    beta = 0.441493616479879;
    
    % Calculating J
    J = A*s1^alpha*exp(-B*s1^beta); % Error less than 0.1% in this range.
    
    % Fit parameters for Delta  
    % Exactly the same as the 12 to 17 range
    C = 0.928133240959474;
    D = 2.326438686265068;
    gamma = 0.720086338581454;
    
    % Calculating Delta
    Delta = 0.5*(s2*(la1/la2)^2)*C*exp(-D*s1^(-gamma)); % Error less than 0.011% in this range.
    
elseif ((s1>=4)&&(s1<8))
    % s1 between 4Er and 8Er
    % Switched to a polynomial model for J
    
    % Fit parameters for J
    p1J = -1.662062243429262e-04;
    p2J = 0.004830606295012;
    p3J = -0.053026148817347;
    p4J = 0.230917549428852;
    
    
    % Calculating J
    J = p1J*s1^3 + p2J*s1^2 + p3J*s1 + p4J; % Error less than 0.1% in this range.
    
    
    % Fit parameters for Delta
    C = 0.961020321737350;
    D = 2.301311224652784;
    gamma = 0.683799690070277;
    
    % Calculating Delta
    Delta = 0.5*(s2*(la1/la2)^2)*C*exp(-D*s1^(-gamma)); % Error less than 0.05% in this range.
    
    
    
elseif ((s1>=1)&&(s1<4))
    % s1 between 1Er and 4Er
    % Switched to polynomial fit for Delta as well.
        
    % Fit parameters for J
    p1J = -2.084465943706971e-04;
    p2J = 0.002510364861444;
    p3J = -0.008029976751638;
    p4J = -0.025758609738462;
    p5J = 0.209662389406680;
    
    % Calculating J
    J = p1J*s1^4 + p2J*s1^3 + p3J*s1^2 + p4J*s1 + p5J; % Error less than 0.05% in this range.
    
    
    % Fit parameters for Delta
    p1D = 9.655792278379040e-05;
    p2D = -1.171970447640949e-04;
    p3D = -0.014162093839709;
    p4D = 0.153946847981079;
    p5D = -0.012252997704736;
    
    % Calculating Delta
    Delta = 0.5*(s2*(la1/la2)^2)*(p1D*s1^4 + p2D*s1^3 + p3D*s1^2 + p4D*s1 + p5D);
    
    
else
    error('The range of the piecewise fit is only determined for s1 between 1Er and 40Er')
end
    
end

