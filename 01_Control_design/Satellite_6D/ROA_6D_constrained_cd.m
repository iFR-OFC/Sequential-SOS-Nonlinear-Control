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
clear
clc
close all

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
um = [-1 -1 -1]'*1.2;
uM = [ 1  1  1]'*1.2;

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
startBuildTime = tic;
% setup first solver
sos = struct('x',[s1;s2;s3;s4],... % decision variables
    'p',[V;kappa;b]);    % parameter

% SOS constraints
sos.('g') = [
    s1;
    s2;
    s3;
    s4;
    s1*(V-b) - nabla(V,x)*subs(f,u,kappa);
    s2*(V-b) + kappa - um;
    s3*(V-b) + uM - kappa;
    s4*(V-b) - g
    ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
opts.error_on_fail = 0; % we do our own check in the main loop

% solver for multiplier
S1 = casos.sossol('S1','mosek',sos,opts);


% setup second solver
cost = dot(g-(V-b),g-(V-b));

sos = struct('x',[V;kappa;b],...            % decision variables
    'f',cost,...
    'p',[s1;s2;s3;s4]);  % parameter

% SOS constraints
sos.('g') = [V-l;
            s1*(V-b) - nabla(V,x)*subs(f,u,kappa);
            s2*(V-b) + kappa - um;
            s3*(V-b) + uM - kappa;
            s4*(V-b) - g
            ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
opts.error_on_fail = 0; % we do our own check in the main loop

% setup solver: solveLyapunov function
S2 = casos.sossol('S2','mosek',sos,opts);
totalBuildTime = toc(startBuildTime);



%% solve problem
fval_old = [];

b = 0.1; % fixed level set

sol2_x = [x'*P0*x;
            -K0*x;
            b];

startVs = tic;
for iter = 1:100

    sol1 = S1('p',sol2_x);

    switch(S1.stats.UNIFIED_RETURN_STATUS)
        case {'SOLVER_RET_SUCCESS'}
            sol1_x = sol1.x;
        otherwise
            fprintf('Solver infeasible in s-step in iteration %d\n',iter)
            break; % in case we are infeasible, leave loop

    end



    sol2 = S2('p',sol1_x );
    switch(S2.stats.UNIFIED_RETURN_STATUS)
        case {'SOLVER_RET_SUCCESS'}
            sol2_x = sol2.x;
        otherwise
            fprintf('Solver infeasible in V step in iteration %d\n',iter)
            break; % in case we are infeasible, leave loop

    end


    % check convergence of cost function
    if ~isempty(fval_old)
        if abs(full(sol2.f-fval_old)) <= 1e-4
            break
        else
            fval_old = sol2.f;
        end
    else
        fval_old = sol2.f;
    end

    % show progress
    fprintf('Iteration %d: f = %g\n',iter,full(sol2.f));

end

totalSolveTime = toc(startVs);
fprintf('----------------------------------- \n',totalBuildTime)
fprintf('Total solve time is %.2f s\n',totalSolveTime)
fprintf('Total build time is %.2f s\n',totalBuildTime)
fprintf('----------------------------------- \n',totalBuildTime)
fprintf('Total time is %.2f s\n',totalBuildTime+totalSolveTime)

%% plotting
figure

% re-scale solution
xd = D*x;

Vfun = to_function(subs(sol2.x(1),x,xd));
gfun = to_function(subs(g,x,xd));

fcontour(@(x2,x3) full(Vfun(0,x2,x3,0) ), [-1 1 -4 4 ], 'b-', 'LevelList', full(sol2.x(end)))
hold on
fcontour(@(x2,x3)  full(gfun(0,x2,x3,0) ), [-1 1 -4 4 ], 'r-', 'LevelList', 0)
hold off
legend('Lyapunov function','Safe set function')
