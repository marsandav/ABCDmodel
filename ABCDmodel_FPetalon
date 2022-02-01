%% Definition of the system with embedded Fabry-Perot etalon

% Illumination beam 
spotsize = 30e-6;                                   % Spot size [m] 
w00 = spotsize/2;                                   % Beam waist [m] 
dlambda = 0.001;                                    % Wavelength step [nm] 
lambda_ini = 1553;                                  % Initial wavelength [nm]
lambda_end = 1555;                                  % Final wavelength [nm]
lambda_vec = (lambda_ini:dlambda:lambda_end)*1e-9;	% Wavelength range [m]
n_0 = 1;                                            % Refractive Index surrounding medium (air)

% FP etalon
n = 1.444;                                          % Refractive Indext spacer (fused silica)
L = 102e-6;                                         % Spacer thickness [m]
R1=0.99;                                            % Reflectivity first mirror
R2=0.99;                                            % Reflectivity second mirror
A1 = 0;                                             % Losses from mirror 1 
A2 = 0;                                             % Losses from mirror 2 
phase_R1f = pi;                                     % Phase shifts from mirror 1 forward

% Spatial sampling of the field 
nradialpts = 1e3;                                   % num points in discretisation across beam 
m_max = 2e3;                                        % maximum iterations
ch_tol = 1e-3;                                      % tolerance to test convergence 
wm_max = 4*L;                                       % maximum extent of integration 

%% ABCD model 

% Initial parameters
rad_vec=linspace(0,1,nradialpts)*wm_max;            % Radius vector
dr=rad_vec(2)-rad_vec(1);                           % Resolution
[rad,lambda] = meshgrid(rad_vec,lambda_vec);        % Radius and wavelength matrices
previous_ITF = ones(size(lambda_vec));              % Set ITF initially to 1

% Mirrors coefficients 
r1=sqrt(R1);                                        % Mirror 1 - Reflection coefficient
t1=sqrt(1-R1-A1);                                   % Mirror 1 - Transmission coefficient
r2=sqrt(R2);                                        % Mirror 2 - Reflection coefficient

% Beam parameters
z00=pi*n_0*w00^2./lambda;                           % Raleigh range
R0=Inf;                                             % Initial radius
k=2*n*pi./lambda;                                   % Wavenumber spacer
kb=2*n_0*pi./lambda;                                % Wavenumber initial medium

% Incident beam
Uin=exp(-rad.^2./w00.^2).*exp(-1i*kb.*rad.^2/2./R0);           

% q-parameter of the incident Gaussian beam
q0=1./(-1i*pi*n*w00^2./lambda./(z00.^2));      

% Initial reflection
Uout=r1*exp(1i*(phase_R1b)).*Uin;  

% 
% ABCD matrices
%
    
% Refraction at a planar boundary (M1)
A=1;                            B=0;	
C=0;                            D=n_0/n;
M_refr=[A B;C D];

% Reverse refraction at a planar boundary (M1)
A=n/n_0;                        B=0;	
C=0;                            D=1;
M_refr_rev=[A B;C D];

% Propagation through homogeneous medium
A=1;                            B=L/n;	
C=0;                            D=1;
M_prop=[A B;C D];

% Reflection from a planar mirror (M1)
A=1;                            B=0;	
C=0;                            D=1;
M_refl1=[A B;C D];

% Reflection from a planar mirror (M2)
A=1;                            B=0;	
C=0;                            D=1;
M_refl2=[A B;C D];

% Round-trip loop
ch = inf;
m=1;
while (m<m_max) && ((m < 100) || (ch > ch_tol))     
          
    M=M_refr_rev*M_prop*M_refl2*M_prop*(M_refl1*M_prop*M_refl2*M_prop)^(m-1)*M_refr;
        
    % Calculate mth field parameters
    A=M(1,1);       B=M(1,2);   
    C=M(2,1);       D=M(2,2);    
    qm=(A*q0+B)./(C*q0+D);                              % mth q-parameter
    qinv=1./qm;                                         % Inverse of q-parameter    
    Rm=1./real(qinv);                                   % Radius of curvature of beam wavefront
    wm=sqrt(lambda./n./pi./(-imag(qinv)));              % Beam width
    w0m=wm./sqrt(1+(pi*n.*wm.^2./lambda./Rm).^2);       % Beam waist
    z0m=pi*n.*w0m.^2./lambda;                           % Rayleigh range
    zm=2*m*L;                                           % Travel path
        
    % Calculate mth U component
    Um=t1*r2*(r1*exp(1i*phase_R1f)*r2)^(m-1)*t1*...
        w00./wm.*exp(-rad.^2./wm.^2).*exp(-1i*k.*zm-1i*k.*rad.^2/2./Rm+1i*atan(n.*pi.*(wm).^2./lambda./Rm)); % mth U component
    
    % Total reflected field complex amplitude
    Uout=Uout+Um; 
    
    % Calculate mth ITF
    Iin=sum(Uin.*conj(Uin).*(2*pi).*rad.*dr,2);
    Iout=sum(Uout.*conj(Uout).*(2*pi).*rad.*dr,2);
    ITF=Iout./Iin;

    % Calculate convergence
    variation = 100*max(abs((ITF-previous_ITF)./previous_ITF));
        
    % Update parameters
    previous_ITF = ITF;    
    ch = variation;    
    m=m+1;
 
end

% Calculate ITF
Iin=sum(Uin.*conj(Uin).*(2*pi).*rad.*dr,2);
Iout=sum(Uout.*conj(Uout).*(2*pi).*rad.*dr,2);
ITF=Iout./Iin;

figure('color','w','DefaultAxesFontSize',14);
plot(lambda_vec*1e9,ITF,'LineWidth',2);shg
ylim([0 1]);
title('ITF');
ylabel('I_r');
xlabel('\lambda (nm)');
legend('Reflection - ABCD model','Location','southwest')
grid on
