%% example1_timing_screening_vs_no_screening.m
% Example 1: timing comparison for the self-transport slice.
% No screening: solve 101 complete H2 synthesis SDPs.
% Screening: solve one recovery-core certificate SDP (48).
% Elapsed times include YALMIP model construction and SDPT3 solving.
% Requirements: MATLAB, YALMIP, SDPT3.

clear; clc;

alpha_grid = logspace(-2,2,101);
epsLMI = 1e-8;  % same numerical margin as the synthesis scripts
opt = sdpsettings('solver','sdpt3','verbose',0);

A = [0.1975 0.7375 0.8957 0.5437 0.7883 0.4944;
     0.8620 0.3033 0.2741 0.3867 0.8529 0.0312;
     0.1256 0.0430 0.9979 0.8222 0.4856 0.8276;
     0.6456 0.8355 0.8344 0.5953 0.8738 0.2643;
     0.4373 0.3693 0.7913 0.7816 0.3338 0.6775;
     0.6126 0.6613 0.6566 0.9877 0.2053 0.7939];

B2 = [0.6818 0.2091 0.0875 0.5936;
      0.6508 0.2725 0.3089 0.0438;
      0.2373 0.7758 0.2309 0.4249;
      0.4774 0.3314 0.9092 0.5216;
      0.9364 0.6031 0.9369 0.8403;
      0.2411 0.1840 0.0319 0.6250];

B1 = eye(6);
C  = [eye(6); zeros(4,6)];
D  = [zeros(6,4); eye(4)];

%% 1) No screening: 101 complete H2 synthesis SDPs
full_status = zeros(numel(alpha_grid),1);
tic;
for k = 1:numel(alpha_grid)
    full_status(k) = solve_self_h2(alpha_grid(k),A,B2,B1,C,D,opt,epsLMI);
end
t_full = toc;

%% 2) Proposed screening: one certificate SDP
tic;

yalmip('clear');
Q = sdpvar(6,6,'symmetric');

N = B2'*Q;
N11 = N(1:2,1:3);
N12 = N(1:2,4:6);
N22 = N(3:4,4:6);

M = A'*Q;
M11 = M(1:3,1:3);
M22 = M(4:6,4:6);

Msum = M11 + M22;
S = (Msum + Msum')/2;

con = [Q >= 0, trace(Q) == 1, N11 + N12 + N22 == 0, Msum == Msum', S >= 0];

sol = optimize(con,[],opt);
t_screen = toc;

%% Report
fprintf('\nNo screening:      101 full H2 SDPs, time = %.4f s\n',t_full);
fprintf('Proposed screening: 1 certificate SDP, time = %.4f s\n',t_screen);
fprintf('Time reduction factor = %.2f x\n',t_full/t_screen);
fprintf('Certificate status: %s\n',yalmiperror(sol.problem));

% Status checks are outside the timed regions.
% YALMIP status 1 means infeasible; status 0 means solved.
statuses = unique(full_status);
for j = 1:numel(statuses)
    fprintf('Full SDP status %d (%s): %d cases\n',statuses(j), ...
        yalmiperror(statuses(j)),nnz(full_status==statuses(j)));
end
if any(full_status~=1)
    warning('Not all complete synthesis SDPs returned the expected infeasible status.');
end
if sol.problem~=0
    warning('Certificate SDP was not reported feasible; check solver output.');
end

%% ------------------------------------------------------------------------
function status = solve_self_h2(alpha,A,B2,B1,C,D,opt,epsLMI)
yalmip('clear');

n = size(A,1);
ell = size(B1,2);
p = size(C,1);
P = sdpvar(n,n,'symmetric');
Z = sdpvar(ell,ell,'symmetric');
X0 = sdpvar(3,3,'full');
R0 = sdpvar(2,3,'full');
gamma = sdpvar(1);

X = [X0 zeros(3); zeros(3) X0];
R = [R0 R0; zeros(2,3) R0];

ZA = A*X-B2*R;
ZC = C*X-D*R;
Qm = [ZA; -X; ZC]*[eye(n),alpha*eye(n),zeros(n,p)];
M0 = [zeros(n),P,zeros(n,p);
      P,zeros(n),zeros(n,p);
      zeros(p,n),zeros(p,n),-eye(p)];
Psi = -M0-Qm-Qm';
L2 = [Z,B1';B1,P];

con = [Psi >= epsLMI*eye(2*n+p), ...
       L2 >= epsLMI*eye(ell+n), ...
       P >= epsLMI*eye(n), Z >= epsLMI*eye(ell), ...
       trace(Z) <= gamma];
sol = optimize(con,gamma,opt);
status = sol.problem;
end
