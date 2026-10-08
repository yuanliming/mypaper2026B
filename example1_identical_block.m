%% example1_identical_block.m
% Example 1: identical-block H2 benchmark
% Reproduces the numerical values reported in the manuscript:
%   - self-transport recovery-core certificate;
%   - non-self transport H2 synthesis over 101 alpha values;
%   - realized H2 costs of the proposed and external controllers.
%
% Requires YALMIP and SDPT3.

clear; clc;

alpha_grid = logspace(-2,2,101);
opt = sdpsettings('solver','sdpt3','verbose',0);

%% Plant and performance data
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

fprintf('Open-loop spectral abscissa = %.4f\n',max(real(eig(A))));

%% 1) Self-transport certificate
% This is the explicit feasible certificate used in the manuscript.
Q = [ 0.048461317856 -0.045182245248  0.075403576679 -0.087259157036  0.076551571556 -0.138534855610;
     -0.045182245248  0.046377851642 -0.064783333968  0.087497825671 -0.078749832468  0.133622179273;
      0.075403576679 -0.064783333968  0.137327835402 -0.124190292918  0.106950848240 -0.207447883314;
     -0.087259157036  0.087497825671 -0.124190292918  0.190227308615 -0.168989204198  0.278193351838;
      0.076551571556 -0.078749832468  0.106950848240 -0.168989204198  0.151399753926 -0.245189807969;
     -0.138534855610  0.133622179273 -0.207447883314  0.278193351838 -0.245189807969  0.426205932560];

Q = (Q+Q')/2;
N = B2'*Q;
M = A'*Q;

Nsum = N(1:2,1:3) + N(1:2,4:6) + N(3:4,4:6);
Msum = M(1:3,1:3) + M(4:6,4:6);
Msum = (Msum+Msum')/2;

fprintf('\nSelf-transport certificate:\n');
fprintf('  max |N11+N12+N22| = %.2e\n',max(abs(Nsum(:))));
fprintf('  lambda_min(Q)      = %.2e\n',min(eig(Q)));
fprintf('  lambda_min(M11+M22)= %.2e\n',min(eig(Msum)));
fprintf('  => self-transport H2 synthesis is infeasible for all alpha > 0.\n');

%% 2) Non-self transport: search the common 101-point alpha grid
best.gamma = inf;

for k = 1:numel(alpha_grid)
    s = solve_nonself(alpha_grid(k),A,B2,B1,C,D,opt);
    if s.ok && s.gamma < best.gamma
        best = s;
        best.alpha = alpha_grid(k);
    end
end

K0 = best.R(1:2,1:3) / best.X(1:3,1:3);

fprintf('\nNon-self transport:\n');
fprintf('  alpha_best = %.4f\n',best.alpha);
fprintf('  gamma_best = %.4f\n',best.gamma);
fprintf('  realized H2^2 = %.4f\n',best.h2sq);
fprintf('  spectral abscissa = %.4f\n',best.sa);
fprintf('  recovered K0 =\n');
disp(K0);

%% 3) External benchmark [14]
K0_ext = [1.7962 2.7236 -1.5147;
          0.1511 0.1290  5.6098];
K_ext = [K0_ext K0_ext; zeros(2,3) K0_ext];

Acl_ext = A-B2*K_ext;
Ccl_ext = C-D*K_ext;
J2_ext = norm(ss(Acl_ext,B1,Ccl_ext,zeros(10,6)),2)^2;

fprintf('External benchmark [14]:\n');
fprintf('  reported certified bound = 163.2\n');
fprintf('  realized H2^2 = %.4f\n',J2_ext);

%% ------------------------------------------------------------------------
function out = solve_nonself(a,A,B2,B1,C,D,opt)
yalmip('clear');

P  = sdpvar(6,6,'symmetric');
Z  = sdpvar(6,6,'symmetric');
X0 = sdpvar(3,3,'full');
R0 = sdpvar(2,3,'full');
g  = sdpvar(1);

% Non-self slice (47)
X = [X0 -X0; zeros(3) X0];
R = [R0 zeros(2,3); zeros(2,3) R0];

Q = [A*X-B2*R; -X; C*X-D*R] * ...
    [eye(6) a*eye(6) zeros(6,10)];

M1 = -[zeros(6) P zeros(6,10);
       P zeros(6) zeros(6,10);
       zeros(10,6) zeros(10,6) -eye(10)] - Q - Q';

M2 = [Z B1'; B1 P];

sol = optimize([M1 >= 0, M2 >= 0, P >= 0, Z >= 0, trace(Z) <= g],g,opt);

out = struct('ok',false,'gamma',inf,'h2sq',nan,'sa',nan,'X',[],'R',[]);
if sol.problem ~= 0, return; end

Xv = value(X);
Rv = value(R);
Kv = Rv/Xv;
Acl = A-B2*Kv;

if rcond(Xv) < 1e-10 || max(real(eig(Acl))) >= 0, return; end

Ccl = C-D*Kv;
h2sq = norm(ss(Acl,B1,Ccl,zeros(10,6)),2)^2;

out = struct('ok',true,'gamma',value(g),'h2sq',h2sq, ...
             'sa',max(real(eig(Acl))),'X',Xv,'R',Rv);
end
