%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction for the longitudinal motion
%                       of the Nasa Generic Transport Model. To increase
%                       the size of the sublevel set we try to minimize the
%                       squared distance to a defined set. Additionally, we
%                       synthesis a linear control law at the same time.
%                       State and control constraints are considered.
%
%   Reference: Modified problem from:
%              Chakraborty, Abhijit and Seiler, Peter and Balas, Gary J.,
%              Nonlinear region of attraction analysis for flight control
%              verification and validation, Control Engineering Practice,
%              2011, doi: 10.1016/j.conengprac.2010.12.001
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%--------------------------------------------------------------------------
close all
clear
clc

import casos.toolboxes.sosopt.cleanpoly


% system states
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);

%% Hubble telescope parameter
J = diag([31046;77217;78754]);

% simple bounds on rates;
omegaMax1 = 0.5*pi/180;
omegaMax2 = 0.2*pi/180;
omegaMax3 = 0.2*pi/180;

x_low =  [-omegaMax1 -omegaMax2 -omegaMax3]';
x_up  =  [ omegaMax1  omegaMax2  omegaMax3]';


% control constraint; assumption is that the box is inside the full
% torque volume. This is roughly estimated visually.
umin = [-1 -1 -1]'*1.2;
umax = [ 1  1  1]'*1.2;

Dx   = diag([1/(x_up(1)-x_low(1)),1/(x_up(2)-x_low(2)),1/(x_up(3)-x_low(3)),1,1,1]); 

Dxin = inv(Dx);

%% dynamics
% cross-product matrix
cpm = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];

% dynamics
B = @(sigma) (1-sigma'*sigma)*eye(3)+ 2*cpm(sigma)+ 2*sigma*sigma';

f =  [-J\cpm(x(1:3))*J*x(1:3) + J\u;
      1/4*B(x(4:6))*x(1:3)]; % omega_dot

% trim point
x0    = [0 0 0 0 0 0]';
u0    = [0,0,0]';

A0 = full(casos.PD(subs(nabla(f,x),[x;u],[x0;u0])));
B0 = full(casos.PD(subs(nabla(f,u),[x;u],[x0;u0])));

% cost function weights
Q = diag([1, 1, 1, 1, 1 ,1]);
R = eye(3)*1;

% generate an initial guess 
[K0,P0] = lqr(full(A0),full(B0),Q,R);

% scaled initial guess for terminal penalty (Lyapunov linear system)
Wval = (inv(Dx)*x)'*P0*(inv(Dx)*x);

% scale dynamics
f = Dx*subs(f,[x;u],[Dx\x;u]);

% state constraint
n = 2;
g0 = (x(1)^2/omegaMax1^2)^(n/2) + (x(2)^2/omegaMax2^2)^(n/2) + (x(3)^2/omegaMax3^2)^(n/2) + ...
     (x(4)^2/0.57^2)^(n/2) + (x(5)^2/0.57^2)^(n/2) + (x(6)^2/0.57^2)^(n/2) - 1;

% re-scale input of state constraints
g = subs(g0,x,Dx\x); 
% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:4));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,4));
s2    = casos.PS.sym('s2',monomials(x,0),[3 1]);
s3    = casos.PS.sym('s3',monomials(x,0),[3 1]);
s4    = casos.PS.sym('s4',monomials(x,0:2));

kappa = casos.PS.sym('k',monomials(x,1),[3 1]);
b     = casos.PS.sym('b');

% enforce positivity
l = 1e-6*(x'*x);


%% setup solver
% options
opts = struct('sossol','mosek');
opts.verbose  = 1;
opts.max_iter = 100;



% cost
cost = dot(g-(V-b),g-(V-b)) ;

sos = struct('x',[V; s1;s2;s3;s4;kappa;b],... % decision variables
    'f',cost, ...                   % cost
    'p',[]);                        % parameter

% SOS constraints
sos.('g') = [s1;
    s2;
    s3;
    s4;
    V-l;
    s1*(V-b)-nabla(V,x)*subs(f,u,kappa)-l;
    s2*(V-b) + kappa - umin;
    s3*(V-b) + umax- kappa;
    s4*(V-b) - g
    ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
opts.sossol_options.newton_solver =[];

% setup solver
S = casos.nlsossol('S1','sequential',sos,opts);

%% solve problem

% initial guess
V0  = g;
s20 = (x'*x)^2;
s30 = [1;1;1];
s40 = [1;1;1];
s50 = x'*x;
K0  = -K0*x;
b0  = 1;

x0 = casos.PD([V0;s20;s30;s40;s50;K0;b0]);

% solve problem
sol = S('x0' ,x0);

% casos.postProcessSolver(S,true);

S.stats

S.stats.single_iterations{end}.conic
