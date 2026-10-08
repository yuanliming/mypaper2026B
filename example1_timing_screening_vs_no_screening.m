%% example1_timing_simple.m
% Example 1: timing comparison
% No screening: solve 101 complete H2 SDPs for the self-transport slice.
% Screening: solve one recovery-core certificate SDP.

clear; clc;

alpha_grid = logspace(-2,2,101);
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
tic;
for k = 1:numel(alpha_grid)
    solve_self_h2(alpha_grid(k),A,B2,B1,C,D,opt);
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
S = M(1:3,1:3) + M(4:6,4:6);

con = [Q >= 0, trace(Q) == 1, N11 + N12 + N22 == 0, S == S', (S + S')/2 >= 0];

sol = optimize(con,[],opt);
t_screen = toc;

%% Report
fprintf('\nNo screening:      101 full H2 SDPs, time = %.4f s\n',t_full);
fprintf('Proposed screening: 1 certificate SDP, time = %.4f s\n',t_screen);
fprintf('Time reduction factor = %.2f x\n',t_full/t_screen);
fprintf('Certificate status: %s\n',yalmiperror(sol.problem));

%% ------------------------------------------------------------------------
function solve_self_h2(a,A,B2,B1,C,D,opt)
yalmip('clear');

P  = sdpvar(6,6,'symmetric');
Z  = sdpvar(6,6,'symmetric');
X0 = sdpvar(3,3,'full');
R0 = sdpvar(2,3,'full');
g  = sdpvar(1);

X = [X0 zeros(3); zeros(3) X0];
R = [R0 R0; zeros(2,3) R0];

Q = [A*X-B2*R; -X; C*X-D*R] * [eye(6) a*eye(6) zeros(6,10)];

M1 = -[zeros(6) P zeros(6,10);
       P zeros(6) zeros(6,10);
       zeros(10,6) zeros(10,6) -eye(10)] - Q - Q';

M2 = [Z B1'; B1 P];

optimize([M1 >= 0, M2 >= 0, P >= 0, Z >= 0, trace(Z) <= g],g,opt);
end
