%% example1_corollary1_robust.m
% Numerically safer Example 1 implementation of Corollary 1.
%
% All modes start from the SAME fixed-coefficient matrix M in (13)-(14).
%   'nullspace'       : exact elimination vec(X,R)=null(M)*theta (recommended)
%   'independentRows' : keep independent full X,R variables, only 45
%                       independent equalities selected from M
%   'fullQR'          : all 126 equalities with YALMIP QR elimination
%
% No hand-coded block restrictions on X or R are used in any mode.
% Requires: MATLAB, YALMIP, SDPT3 
clear; clc;
mode = 'nullspace';  % alternatively 'independentRows' or 'fullQR'
alpha = 10^(-0.04);
epsLMI = 1e-8;

%% Original identical-block plant and H2 output
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
n = 6; m = 4; ell = 6; p = n+m;
B1 = eye(n);
C = [eye(n); zeros(m,n)];
D = [zeros(n,m); eye(m)];

%% Controller direction space and fixed non-self reference
J = [1 1; 0 1];
Kbasis = zeros(m*n,6);
t = 0;
for c = 1:3
    for r = 1:2
        t = t+1;
        E = zeros(2,3);
        E(r,c) = 1;
        N = kron(J,E);
        Kbasis(:,t) = N(:);
    end
end
H = null(Kbasis')';
h = zeros(size(H,1),1);
F = kron(J\eye(2),eye(3));
[M,GF] = build_corollary1_encoding(H,h,F,m);
rM = rank(M);
V = null(M);
fprintf('Corollary 1: M=%d x %d, rank=%d, nullity=%d\n', ...
        size(M,1),size(M,2),rM,size(V,2));
assert(rM==45 && size(V,2)==15,'Unexpected factor-space dimension.');
assert(norm(M*V,inf)<1e-10,'Numerically unstable null space of M.');

%% Convex variables: choose one equivalent representation of ker(M)
yalmip('clear');
P = sdpvar(n,n,'symmetric');
Z = sdpvar(ell,ell,'symmetric');
gamma = sdpvar(1);

switch mode
    case 'nullspace'
        % Corollary 1 Eq. (14), analytically eliminating the equalities.
        theta = sdpvar(size(V,2),1);
        v = V*theta;
        X = reshape(v(1:n*n),n,n);
        R = reshape(v(n*n+1:end),m,n);
        factor_eq = [];
        opts = sdpsettings('solver','sdpt3','verbose',1);
    case 'independentRows'
        % Still full and independent matrix variables X,R, but no
        % redundant copies of the 45 equality constraints.
        X = sdpvar(n,n,'full');
        R = sdpvar(m,n,'full');
        [~,~,piv] = qr(M','vector');
        Mred = M(piv(1:rM),:);
        assert(rank(Mred)==rM);
        assert(norm(M*null(Mred),inf)<1e-10);
        factor_eq = (Mred*[X(:);R(:)] == 0);
        opts = sdpsettings('solver','sdpt3','verbose',1, ...
                           'removeequalities',1);
    case 'fullQR'
        % Same full 126 equalities as the original script; ask YALMIP
        % to remove them by QR before sending the SDP to SDPT3.
        X = sdpvar(n,n,'full');
        R = sdpvar(m,n,'full');
        factor_eq = (M*[X(:);R(:)] == 0);
        opts = sdpsettings('solver','sdpt3','verbose',1, ...
                           'removeequalities',1);
    otherwise
        error('Unknown mode: %s',mode);
end

%% Exactly the same H2 LMI as in the original example
ZA = A*X-B2*R;
ZC = C*X-D*R;
Qm = [ZA;-X;ZC]*[eye(n),alpha*eye(n),zeros(n,p)];
M0 = [zeros(n),P,zeros(n,p);
      P,zeros(n),zeros(n,p);
      zeros(p,n),zeros(p,n),-eye(p)];
Psi = -M0-Qm-Qm';
L2 = [Z,B1';B1,P];
con = [factor_eq, ...
       Psi >= epsLMI*eye(2*n+p), ...
       L2 >= epsLMI*eye(ell+n), ...
       P >= epsLMI*eye(n), ...
       Z >= epsLMI*eye(ell), ...
       trace(Z) <= gamma];

sol = optimize(con,gamma,opts);
fprintf('method=%s, problem=%d, info=%s\n',mode,sol.problem,sol.info);
if sol.problem~=0
    warning('SDP did not return an optimum. Inspect SDPT3 log.');
    return;
end

%% Independent residual and closed-loop validation
Xv = value(X); Rv = value(R);
Pv = value(P); Zv = value(Z);
Psiv = value(Psi); L2v = value(L2);
if any(~isfinite(Xv(:))) || rcond(Xv)<1e-10
    error('The recovered denominator X is invalid or ill-conditioned.');
end
Kv = Rv/Xv;
Acl = A-B2*Kv;
assert(max(real(eig(Acl)))<0,'Recovered controller is not Hurwitz.');
J2 = norm(ss(Acl,B1,C-D*Kv,zeros(p,ell)),2)^2;
eqres = norm(M*[Xv(:);Rv(:)],inf);
Kres = norm(H*Kv(:)-h,inf);
fprintf('alpha=%.7f  gamma=%.7f  H2^2=%.7f\n',alpha,value(gamma),J2);
fprintf('factor residual %.3e, gain equality residual %.3e\n',eqres,Kres);
fprintf('rcond(X)=%.3e  spectral abscissa(Acl)=%.6f\n', ...
        rcond(Xv),max(real(eig(Acl))));
fprintf('min eig Psi=%.3e, L2=%.3e, P=%.3e, Z=%.3e\n', ...
        min(eig((Psiv+Psiv')/2)), ...
        min(eig((L2v+L2v')/2)), ...
        min(eig((Pv+Pv')/2)), ...
        min(eig((Zv+Zv')/2)));
fprintf('H2 bound gap = %.6f\n',value(gamma)-J2);
disp('Recovered K0 (top-left 2 x 3 block):');
disp(Kv(1:2,1:3));

function [M,G,Nbar] = build_corollary1_encoding(H,h,F,m)
%BUILD_COROLLARY1_ENCODING Construct the fixed matrix in Corollary 1.
%   [M,G,Nbar] = build_corollary1_encoding(H,h,F,m)
%   encodes, for X invertible,
%     S_H * X = S_G   and   H*vec(R/X) = h
%   as the fixed-coefficient linear constraint
%     M * [X(:); R(:)] == 0.
%
%   H: q-by-(m*n), full-row-rank controller equality matrix
%   h: q-by-1 affine offset
%   F: n-by-n fixed nonsingular slice reference, S_G = S_H*F
%   m: number of controller rows
%
%   G = H*kron(F^(-T),eye(m)), as in Eq. (20).
%   M is assembled exactly as Eq. (13), using null([H,-h]).
%   Neither M nor G contains optimization variables.
%
%   Put this file on the MATLAB path along with the example scripts.

    [q,mn] = size(H);
    n = size(F,1);
    assert(isequal(size(F),[n,n]),'F must be square.');
    assert(mn == m*n,'size(H,2) must equal m*n.');
    assert(numel(h)==q,'h has incorrect dimension.');
    h = h(:);
    assert(rank(H)==q,'H must have full row rank.');
    assert(rcond(F)>1e-12,'F must be nonsingular.');

    % Eq. (20): use a linear solve rather than explicitly forming inv(F').
    G = H*kron(F'\eye(n),eye(m));

    % Eqs. (11)-(12): columns are [vec(N_j); eta_j].
    Nbar = null([H,-h]);
    dplus1 = mn+1-q;
    assert(size(Nbar,2)==dplus1,'Unexpected null-space dimension.');

    % Eq. (13): M is [q*(d+1)]-by-[n*n+m*n].
    M = zeros(q*dplus1,n*n+mn);
    for j=1:dplus1
        Nj = reshape(Nbar(1:mn,j),m,n);
        eta = Nbar(end,j);
        rows = (j-1)*q+(1:q);
        M(rows,1:n*n) = G*kron(eye(n),Nj);
        M(rows,n*n+1:end) = -eta*G;
    end
end
