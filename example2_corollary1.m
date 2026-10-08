%% example2_corollary1.m
% Mixed-affine H2 synthesis using the FIXED-COEFFICIENT EQUALITY (14).
% Replaces the xi,c1,c2 parameterization (54) by independent YALMIP
% decision matrices X (3x3), R (2x3), coupled by M*[X(:);R(:)] == 0.
% Same plant, references, alpha grid, performance LMI and objective as
% the original example2_mixed_affine.m (paper Eqs. (41)-(42)).
%
% Requirements: MATLAB, YALMIP, SDPT3
%
% Set RUN_FULL_GRID=false for a quick single-case demonstration.
% Set RUN_FULL_GRID=true to evaluate all 10 references x 101 alpha values.

clear; clc;
RUN_FULL_GRID = true;
alpha_grid = logspace(-2,2,101);
opt = sdpsettings('solver','sdpt3','verbose',0);

%% Same data as Eqs. (50)-(53)
A  = [-1 2 1; 2 -1 0; 4 4 -2];
B2 = [2 0; 0 2; 2 2];
B1 = eye(3);
C  = [eye(3); zeros(2,3)];
D  = [zeros(3,2); eye(2)];

% vec(K) = [k11 k21 k12 k22 k13 k23]'.
H = [1 0 0 0 0 0;
     0 0 0 0 0 1;
     0 -1 1 0 0 0;
     0 0 0 1 2 0];
h = [0; 1; 0; 3];

% Only used as a numerical witness to TEST M, not to CONSTRUCT M.
Kp = [0 0 0; 0 3 1];
N1 = [0 1 0; 1 0 0];
N2 = [0 0 1; 0 -2 0];
assert(norm(H*Kp(:)-h,inf) < 1e-12);

%% Same reference library as Eq. (55)
Fb = diag([1 1 2]);
E12 = [0 1 0; 0 0 0; 0 0 0];
E23 = [0 0 0; 0 0 1; 0 0 0];
F = cell(10,1); label = cell(10,1); tau = nan(10,1);
F{1}=eye(3); F{2}=Fb;
label{1}='F1'; label{2}='F2';
k = 2;
for t = [0.2 0.5]
    F{k+1} = Fb*(eye(3)+t*E12); label{k+1}='F12+'; tau(k+1)=t;
    F{k+2} = Fb*(eye(3)-t*E12); label{k+2}='F12-'; tau(k+2)=t;
    F{k+3} = Fb*(eye(3)+t*E23); label{k+3}='F23+'; tau(k+3)=t;
    F{k+4} = Fb*(eye(3)-t*E23); label{k+4}='F23-'; tau(k+4)=t;
    k=k+4;
end

% F10 is the best-bound reference reported in Table III.
if RUN_FULL_GRID
    ref_ids = 1:10;
    a_values = alpha_grid;
else
    ref_ids = 10;
    a_values = 10^(-0.64);  % approx. 0.2291
end

best = repmat(struct('ok',false,'gamma',inf,'h2',nan, ...
              'alpha',nan,'K',[],'eqres',nan),10,1);
gamma_all = nan(10,numel(a_values));
fprintf('Corollary 1: full factor-space equality formulation\n');

for j=ref_ids
    [M,GF] = build_corollary1_encoding(H,h,F{j},2);
    assert(size(M,2)==15);
    assert(rank(M)==12,'Expected 3-dimensional factor space.');

    % Explicit algebraic consistency checks (no optimization here).
    test_factors = [1.2*Kp - 0.7*N1 + 0.4*N2]*F{j};
    Xtest = 1.2*F{j};
    assert(norm(M*[Xtest(:);test_factors(:)],inf)<1e-10);
    assert(norm(GF*kron(F{j}',eye(2))*N1(:),inf)<1e-10);

    for ia=1:numel(a_values)
        alpha=a_values(ia);
        s=solve_with_corollary1(alpha,M,A,B2,B1,C,D,H,h,opt);
        if s.ok
            gamma_all(j,ia)=s.gamma;
            if s.gamma<best(j).gamma
                best(j)=s;
                best(j).alpha=alpha;
            end
        end
    end
    fprintf('%-5s tau=%4.1f: alpha=%8.4f gamma=%10.4f H2^2=%9.4f, eqres=%.2e\n', ...
         label{j},tau(j),best(j).alpha,best(j).gamma,best(j).h2,best(j).eqres);
end

if RUN_FULL_GRID
    % The six rows correspond to Table III, selecting the best tau in
    % each signed-shear family.
    groups = {1,2,[3 7],[4 8],[5 9],[6 10]};
    fprintf('\nTable III family winners (Corollary 1 encoding):\n');
    for g=1:numel(groups)
        ids=groups{g};
        [~,idloc]=min(arrayfun(@(i)best(i).gamma,ids));
        j=ids(idloc);
        fprintf('%-5s tau=%4.1f alpha=%.4f gamma=%.4f H2^2=%.4f\n', ...
            label{j},tau(j),best(j).alpha,best(j).gamma,best(j).h2);
    end
end

%% Local solver with direct fixed-coefficient equality
function out=solve_with_corollary1(a,M,A,B2,B1,C,D,H,h,opt)
    yalmip('clear');
    n=3; m=2; ell=3; p=5;
    P=sdpvar(n,n,'symmetric');
    Z=sdpvar(ell,ell,'symmetric');
    X=sdpvar(n,n,'full');
    R=sdpvar(m,n,'full');
    gamma=sdpvar(1);

    % ***** Corollary 1, Eq. (14): FIXED coefficients! *****
    factor_eq = M*[X(:);R(:)] == 0;

    ZA=A*X-B2*R;
    ZC=C*X-D*R;
    Qm=[ZA; -X; ZC]*[eye(n),a*eye(n),zeros(n,p)];
    M0=[zeros(n),P,zeros(n,p);
        P,zeros(n),zeros(n,p);
        zeros(p,n),zeros(p,n),-eye(p)];
    Psi=-M0-Qm-Qm';
    L2=[Z,B1';B1,P];
    epsLMI=1e-8;
    con=[factor_eq, ...
         Psi>=epsLMI*eye(2*n+p), ...
         L2>=epsLMI*eye(ell+n), ...
         P>=epsLMI*eye(n), Z>=epsLMI*eye(ell), ...
         trace(Z)<=gamma];
    sol=optimize(con,gamma,opt);

    out=struct('ok',false,'gamma',inf,'h2',nan, ...
               'alpha',nan,'K',[],'eqres',nan);
    if sol.problem~=0, return; end
    Xv=value(X); Rv=value(R);
    if rcond(Xv)<1e-10, return; end
    Kv=Rv/Xv;
    Acl=A-B2*Kv;
    if max(real(eig(Acl)))>=-1e-8, return; end

    eqres=norm(M*[Xv(:);Rv(:)],inf);
    Kres=norm(H*Kv(:)-h,inf);
    if eqres>1e-5 || Kres>1e-5
        warning('Unacceptably large equality residual: M=%.3e, K=%.3e',eqres,Kres);
        return;
    end
    sys=ss(Acl,B1,C-D*Kv,zeros(p,ell));
    out=struct('ok',true,'gamma',value(gamma), ...
               'h2',norm(sys,2)^2,'alpha',a,'K',Kv,'eqres',eqres);
end

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
