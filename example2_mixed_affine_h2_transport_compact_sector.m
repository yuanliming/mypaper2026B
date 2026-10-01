clear; clc; close all; tic;

%% Mixed affine-equality structured H2 synthesis
% K = [k11 k12 k13; k21 k22 k23]
% Constraints: k11=0, k23=1, k12=k21, k22+2*k13=3.
% vec(K) = [k11;k21;k12;k22;k13;k23].

alpha_grid = logspace(-2,2,101);
opt = sdpsettings('solver','sdpt3','verbose',0);

%% Plant
A  = [-1 2 1; 2 -1 0; 4 4 -2];
B2 = [2 0; 0 2; 2 2];
B1 = eye(3);
C  = [eye(3); zeros(2,3)];
D  = [zeros(3,2); eye(2)];
[n,m] = deal(size(A,1),size(B2,2));

fprintf('Open-loop spectral abscissa = %.6g\n',max(real(eig(A))));

%% Affine controller constraints H*vec(K)=h
H = [1  0  0  0  0  0;
     0  0  0  0  0  1;
     0 -1  1  0  0  0;
     0  0  0  1  2  0];
h = [0;1;0;3];
assert(rank(H)==size(H,1),'H must have full row rank.');

%% Fixed reference library
% F1: self transport; F2: balanced D; F3-F6: signed shears about D.
tau = 0.5;
Dbal = diag([1 1 2]);
Fcell = {eye(3); Dbal; ...
         Dbal*(eye(3)+tau*[0 1 0;0 0 0;0 0 0]); ...
         Dbal*(eye(3)-tau*[0 1 0;0 0 0;0 0 0]); ...
         Dbal*(eye(3)+tau*[0 0 0;0 0 1;0 0 0]); ...
         Dbal*(eye(3)-tau*[0 0 0;0 0 1;0 0 0])};
Fname = {'Self transport','Transport F2','Transport F3', ...
         'Transport F4','Transport F5','Transport F6'};
nF = numel(Fcell);

%% Precompute transported constraints
Gcell = cell(nF,1); Mcell = cell(nF,1); sector_ok = false(nF,1);
for j = 1:nF
    F = Fcell{j};
    if rcond(F)<1e-12, error('Reference F{%d} is singular.',j); end
    % Here A_H = span{I_3}. Thus the denominator-sector test
    % Y in A_H, He(YF) >= I is feasible iff He(F) is definite.
    ev = eig((F+F')/2);
    sector_ok(j) = all(ev>1e-10) || all(ev<-1e-10);
    fprintf('%-18s denominator-sector pass = %d\n',Fname{j},sector_ok(j));
    if ~sector_ok(j), continue; end
    FinvT = (F\eye(n)).';
    Gcell{j} = H*kron(FinvT,eye(m));
    Mcell{j} = affine_transport_encoder(H,h,Gcell{j},m,n);
end

%% Solve all slices on the same alpha grid
gamma_all = nan(nF,numel(alpha_grid));
h2sq_all  = nan(nF,numel(alpha_grid));
best = repmat(struct('ok',false,'gamma',inf,'h2sq',nan,'alpha',nan, ...
    'X',[],'R',[],'K',[],'residual',nan,'margin',nan),nF,1);

for j = 1:nF
    if ~sector_ok(j), continue; end
    fprintf('\nSolving %s ...\n',Fname{j});
    for k = 1:numel(alpha_grid)
        a = alpha_grid(k);
        s = solve_affine_slice(a,Mcell{j},H,h,A,B2,B1,C,D,opt);
        if s.ok
            gamma_all(j,k) = s.gamma;
            h2sq_all(j,k) = s.h2sq;
            if s.gamma < best(j).gamma
                best(j) = s; best(j).alpha = a;
            end
        end
    end
end

%% Summary
fprintf('\n============================================================\n');
fprintf('Search alpha in [%.2e, %.2e], %d points\n', ...
    alpha_grid(1),alpha_grid(end),numel(alpha_grid));
fprintf('============================================================\n');
for j = 1:nF
    if best(j).ok
        fprintf(['%-18s alpha*=%.6g, gamma*=%.6g, actual H2^2=%.6g, ' ...
            'gamma/H2^2=%.6g, residual=%.3e, margin=%.6g\n'], ...
            Fname{j},best(j).alpha,best(j).gamma,best(j).h2sq, ...
            best(j).gamma/best(j).h2sq,best(j).residual,best(j).margin);
    else
        fprintf('%-18s no feasible grid point found.\n',Fname{j});
    end
end

%% Best non-self slice
gammas_nonself = arrayfun(@(s)s.gamma,best(2:end));
[best_gamma_nonself,idx0] = min(gammas_nonself); idx_tr = idx0+1;
fprintf('\n================ BEST NON-SELF SLICE ================\n');
if isfinite(best_gamma_nonself)
    s = best(idx_tr); K = s.K;
    fprintf('Reference       : %s\n',Fname{idx_tr});
    fprintf('alpha*          : %.8g\n',s.alpha);
    fprintf('gamma*          : %.8g\n',s.gamma);
    fprintf('actual H2^2     : %.8g\n',s.h2sq);
    fprintf('gamma/H2^2      : %.8g\n',s.gamma/s.h2sq);
    fprintf('spectral margin : %.8g\n',s.margin);
    fprintf('affine residual : %.3e\n',s.residual);
    disp('F ='); disp(Fcell{idx_tr});
    disp('K ='); disp(K);
    fprintf('k11               = %.10f  (target 0)\n',K(1,1));
    fprintf('k23               = %.10f  (target 1)\n',K(2,3));
    fprintf('k12-k21           = %.3e   (target 0)\n',K(1,2)-K(2,1));
    fprintf('k22+2*k13         = %.10f  (target 3)\n',K(2,2)+2*K(1,3));
else
    fprintf('No non-self transport slice was feasible on the grid.\n');
end

%% Self-transport baseline
fprintf('\n================ SELF-TRANSPORT BASELINE ================\n');
if best(1).ok
    s0 = best(1);
    fprintf('alpha*          : %.8g\n',s0.alpha);
    fprintf('gamma*          : %.8g\n',s0.gamma);
    fprintf('actual H2^2     : %.8g\n',s0.h2sq);
    fprintf('gamma/H2^2      : %.8g\n',s0.gamma/s0.h2sq);
    fprintf('spectral margin : %.8g\n',s0.margin);
    fprintf('affine residual : %.3e\n',s0.residual);
    disp('K_self ='); disp(s0.K);
else
    fprintf('No feasible self-transport point found on the grid.\n');
end

%% Plot certified H2^2 upper bounds
figure; hold on;
for j = 1:nF
    id = ~isnan(gamma_all(j,:));
    if any(id)
        semilogx(alpha_grid(id),gamma_all(j,id),'-o', ...
            'LineWidth',1,'MarkerSize',3,'DisplayName',Fname{j});
    end
end
grid on; xlabel('\alpha'); ylabel('\gamma');
legend('show','Location','best'); title('Certified squared H_2 upper bounds');
toc;

%% Local functions
function M = affine_transport_encoder(H,h,G,m,n)
Hbar = [H -h]; Nbar = null(Hbar);
q = size(H,1); dbar = size(Nbar,2);
M = zeros(q*dbar,n*n+m*n);
for j = 1:dbar
    Nj = reshape(Nbar(1:m*n,j),m,n);
    eta = Nbar(end,j); rows = (j-1)*q+(1:q);
    M(rows,:) = [G*kron(eye(n),Nj), -eta*G];
end
end

function out = solve_affine_slice(a,Menc,H,h,A,B2,B1,C,D,opt)
yalmip('clear');
n = size(A,1); m = size(B2,2); nw = size(B1,2); nz = size(C,1);
P = sdpvar(n,n,'symmetric'); Z = sdpvar(nw,nw,'symmetric');
X = sdpvar(n,n,'full'); R = sdpvar(m,n,'full'); g = sdpvar(1);

Q = [A*X-B2*R; -X; C*X-D*R]*[eye(n),a*eye(n),zeros(n,nz)];
M0 = [zeros(n) P zeros(n,nz); P zeros(n) zeros(n,nz); ...
      zeros(nz,n) zeros(nz,n) -eye(nz)];
L1 = -M0-Q-Q';
L2 = [Z B1'; B1 P];
eq = Menc*[X(:);R(:)]==0; eps_lmi = 1e-8;
con = [eq, L1>=eps_lmi*eye(size(L1)), L2>=eps_lmi*eye(size(L2)), ...
       P>=eps_lmi*eye(n), Z>=eps_lmi*eye(nw), trace(Z)<=g];
sol = optimize(con,g,opt);
out = struct('ok',false,'gamma',inf,'h2sq',nan,'alpha',nan, ...
    'X',[],'R',[],'K',[],'residual',nan,'margin',nan);
if sol.problem~=0, return; end

Xv = value(X); Rv = value(R); gv = value(g);
if ~isfinite(gv) || rcond(Xv)<1e-10, return; end
Kv = Rv/Xv; Acl = A-B2*Kv; sa = max(real(eig(Acl)));
if sa>=-1e-8, return; end
residual = norm(H*Kv(:)-h,inf);
if residual>1e-5, return; end
Ccl = C-D*Kv; sys = ss(Acl,B1,Ccl,zeros(nz,nw));
h2sq = norm(sys,2)^2;
if ~isfinite(h2sq), return; end
out = struct('ok',true,'gamma',gv,'h2sq',h2sq,'alpha',nan, ...
    'X',Xv,'R',Rv,'K',Kv,'residual',residual,'margin',-sa);
end

%%% the optimized gains

best(1).K
best(2).K
best(3).K
best(4).K
best(5).K
best(6).K