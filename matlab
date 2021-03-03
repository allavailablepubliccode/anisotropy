clear;clc;close all;

% specify states and parameters
%==========================================================================
N = 100;                        % number of time points

% model states, phi (x.p), and theta (x.t)
%--------------------------------------------------------------------------
x.p             =  0.1;         % initialise phi at different values
x.p_xp1         =  0.2;         
x.p_xm1         =  0.3;
x.p_yp1         =  0.4;
x.p_ym1         =  0.5;
x.p_xp1_yp1     =  0.6;
x.p_xp1_ym1     =  0.7;
x.p_xm1_yp1     =  0.8;
x.p_xm1_ym1     =  0.9;

x.t             =  0.1;         % initialise theta at different values
x.t_xp1         =  0.2;
x.t_xm1         =  0.3;
x.t_yp1         =  0.4;
x.t_ym1         =  0.5;
x.t_xp1_yp1     =  0.6;
x.t_xp1_ym1     =  0.7;
x.t_xm1_yp1     =  0.8;
x.t_xm1_ym1     =  0.9;

% model parameters - delta = 0 (wave equation)
%--------------------------------------------------------------------------
P.del           = 0; 

% observation function (phi - to generate timeseries)
%--------------------------------------------------------------------------
g = @(x,v,P) normalize([x.p;x.p_xp1;x.p_xm1;x.p_yp1;x.p_ym1;x.p_xp1_yp1;...
    x.p_xp1_ym1;x.p_xm1_yp1;x.p_xm1_ym1]);

% equations of motion
%--------------------------------------------------------------------------
f = @(x,v,P) [normalize([x.t;x.t_xp1;x.t_xm1;x.t_yp1;x.t_ym1;...
    x.t_xp1_yp1;x.t_xp1_ym1;x.t_xm1_yp1;x.t_xm1_ym1]);...
    normalize([(x.p_xp1-2*x.p+x.p_xm1)+x.p.^(2*P.del).*((x.p_yp1-2*x.p+...
    x.p_ym1)+P.del.*0.5*(x.p_yp1-x.p_ym1).^2./x.p);...
    (x.p_xp1-2*x.p+x.p_xm1)+x.p_xp1.^(2*P.del).*((x.p_xp1_yp1-2*...
    x.p_xp1+x.p_xp1_ym1)+P.del.*0.5*(x.p_xp1_yp1-x.p_xp1_ym1).^2./...
    x.p_xp1);(x.p_xm1-2*x.p+x.p_xp1)+x.p_xm1.^(2*P.del).*((x.p_xm1_yp1-...
    2*x.p_xm1+x.p_xm1_ym1)+P.del.*0.5*(x.p_xm1_yp1-x.p_xm1_ym1).^2./...
    x.p_xm1);(x.p_xp1_yp1-2*x.p_yp1+x.p_xm1_yp1)+x.p_yp1.^(2*P.del).*...
    ((x.p_xp1_yp1-2*x.p_xp1+x.p_xp1_ym1)+P.del.*(x.p_yp1-x.p).^2./...
    x.p_yp1);(x.p_xp1_ym1-2*x.p_ym1+x.p_xm1_ym1)+x.p_ym1.^(2*P.del).*...
    ((x.p_xp1_ym1-2*x.p_xp1+x.p_xp1_yp1)+P.del.*(x.p-x.p_ym1).^2./...
    x.p_ym1);(x.p_xp1_yp1-2*x.p_yp1+x.p_xm1_yp1)+x.p_xp1_yp1.^...
    (2*P.del).*((x.p_xp1_yp1-2*x.p_xp1+x.p_xp1_ym1)+P.del.*(x.p_xp1_yp1-...
    x.p_xp1).^2./x.p_xp1_yp1);(x.p_xp1_ym1-2*x.p_ym1+x.p_xm1_ym1)+...
    x.p_xp1_ym1.^(2*P.del).*((x.p_xp1_ym1-2*x.p_xp1+x.p_xp1_yp1)+...
    P.del.*(x.p_xp1-x.p_xp1_ym1).^2./x.p_xp1_ym1);(x.p_xm1_yp1-2*...
    x.p_yp1+x.p_xp1_yp1)+x.p_xm1_yp1.^(2*P.del).*((x.p_xm1_yp1-...
    2*x.p_xm1+x.p_xm1_ym1)+P.del.*(x.p_xm1_yp1-x.p_xm1).^2./...
    x.p_xm1_yp1);(x.p_xm1_ym1-2*x.p_ym1+x.p_xp1_ym1)+x.p_xm1_ym1.^...
    (2*P.del).*((x.p_xm1_ym1-2*x.p_xm1+x.p_xm1_yp1)+P.del.*(x.p_xm1-...
    x.p_xm1_ym1).^2./x.p_xm1_ym1)])];

% causes or exogenous input (zero perturbation)
%--------------------------------------------------------------------------
U               = zeros(1,N);

% first level state space model
%--------------------------------------------------------------------------
M(1).x          = x;            % initial states 
M(1).f          = f;            % equations of motion
M(1).g          = g;            % observation mapping
M(1).pE         = P;            % model parameters
M(1).V          = exp(26);      % precision of observation noise
M(1).W          = exp(26);      % precision of state noise

% second level causes or exogenous forcing term
%--------------------------------------------------------------------------
M(2).v          = 0;            % initial causes
M(2).V          = exp(26);      % precision of exogenous causes

% create isotropic data with known parameters (P)
%==========================================================================
DEM_iso     	= spm_DEM_generate(M,U,P);

% change delta to non-zero value for anisotropic case
%--------------------------------------------------------------------------
P.del           = 6;
M(1).pE         = P;

% create anisotropic data with known parameters (P)
%==========================================================================
DEM_aniso    	= spm_DEM_generate(M,U,P);

% Now try to recover model parameters from data features
%==========================================================================

% initialization of priors over parameters
%--------------------------------------------------------------------------
DEM_aniso.M(1).pE.del   = 0;        % set prior parameter delta to zero

pC.del                  = 1/64;                 % prior variance 
DEM_iso.M(1).pC         = diag(spm_vec(pC));
DEM_aniso.M(1).pC       = diag(spm_vec(pC));

% Inversion using generalised filtering 
%==========================================================================
LAP_iso             = spm_DEM(DEM_iso);
LAP_ani             = spm_DEM(DEM_aniso);

% use Bayesian model reduction to test different hypotheses
%==========================================================================
model{1} = 'isotropic';
model{2} = 'anisotropic';

% apply precise shrinkage priors to of diagonal coupling elements
%--------------------------------------------------------------------------
PC{1} = pC; PC{1}.del = 0;
PC{2} = pC;

%  evaluate the evidence for these new models or prior constraints
%--------------------------------------------------------------------------
qE    = LAP_iso.qP.P{1};            % isotropic case
qC    = LAP_iso.qP.C;
pE    = LAP_iso.M(1).pE;
pC    = LAP_iso.M(1).pC;
for m = 1:numel(PC)
    rC          = diag(spm_vec(PC{m}));
    F_iso(m)    = spm_log_evidence(qE,qC,pE,pC,pE,rC);
end

qE    = LAP_ani.qP.P{1};                % anisotropic case
qC    = LAP_ani.qP.C;
pE    = LAP_ani.M(1).pE;
pC    = LAP_ani.M(1).pC;
for m = 1:numel(PC)
    rC          = diag(spm_vec(PC{m}));
    F_aniso(m)  = spm_log_evidence(qE,qC,pE,pC,pE,rC);
end

% report marginal log likelihood or evidence
%--------------------------------------------------------------------------
F_iso   = F_iso - min(F_iso);           % isotropic case
F_aniso = F_aniso - min(F_aniso);       % anisotropic case

% plot result of Bayesian model reduction
%--------------------------------------------------------------------------
figure;

subplot(2,2,1);
bar(F_iso,'c')
title('Log evidence - iso','FontSize',16)
xlabel(model), axis square, box off

subplot(2,2,2);
bar(spm_softmax(F_iso(:)),'c')
title('Probability - iso','FontSize',16)
xlabel(model), axis square, box off

subplot(2,2,3);
bar(F_aniso,'c')
title('Log evidence - aniso','FontSize',16)
xlabel(model), axis square, box off

subplot(2,2,4);
bar(spm_softmax(F_aniso(:)),'c')
title('Probability - aniso','FontSize',16)
xlabel(model), axis square, box off

% plot synthetic data
%--------------------------------------------------------------------------
ts = 1.05:1.05:N*1.05; 
x = repmat(cos(ts),[9,1]); 
y = repmat(sin(ts),[9,1]);
figure;
subplot(2,2,1);
plot(LAP_ani.Y'.*x',LAP_ani.Y'.*y');
title('anisotropic')
axis equal;
subplot(2,2,2);
plot(LAP_iso.Y'.*x',LAP_iso.Y'.*y');
axis equal;
title('isotropic')
