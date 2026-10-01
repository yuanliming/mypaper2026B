clear; clc; tic;
alpha_grid = logspace(-2,2,101); %logspacelinspace
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
gamma_fdi = nan(size(alpha_grid));
gamma_tr  = nan(size(alpha_grid));
h2sq_fdi  = nan(size(alpha_grid));
h2sq_tr   = nan(size(alpha_grid));
best_fdi = struct('gamma',inf);
best_tr  = struct('gamma',inf);
for k = 1:numel(alpha_grid)
    a = alpha_grid(k);
    s = solve_slice(a,0,A,B2,B1,C,D,opt);
    if s.ok
        gamma_fdi(k) = s.gamma;
        h2sq_fdi(k) = s.h2sq;
        if s.gamma < best_fdi.gamma, best_fdi = s; best_fdi.alpha = a; end
    end
    s = solve_slice(a,1,A,B2,B1,C,D,opt);
    if s.ok
        gamma_tr(k) = s.gamma;
        h2sq_tr(k) = s.h2sq;
        if s.gamma < best_tr.gamma, best_tr = s; best_tr.alpha = a; end
    end
end
fprintf('\nSearch alpha in [%.2e, %.2e], %d points\n', ...
    alpha_grid(1),alpha_grid(end),numel(alpha_grid));
if isfinite(best_fdi.gamma)
    fprintf('FDI: alpha*=%.6g, gamma*=%.6g, actual H2^2=%.6g, gamma/H2^2=%.6g\n', ...
        best_fdi.alpha,best_fdi.gamma,best_fdi.h2sq,best_fdi.gamma/best_fdi.h2sq);
else
    fprintf('FDI: no feasible point found on the grid.\n');
end
if isfinite(best_tr.gamma)
    fprintf('TR : alpha*=%.6g, gamma*=%.6g, actual H2^2=%.6g, gamma/H2^2=%.6g\n', ...
        best_tr.alpha,best_tr.gamma,best_tr.h2sq,best_tr.gamma/best_tr.h2sq);
    X0 = best_tr.X(1:3,1:3);
    R0 = best_tr.R(1:2,1:3);
    K0 = R0/X0;
    disp('X0 ='); disp(X0);
    disp('R0 ='); disp(R0);
    disp('K0 ='); disp(K0);
    fprintf('spectral abscissa = %.6g\n',max(real(eig(A-B2*best_tr.K))));
end
i1 = ~isnan(gamma_fdi);
i2 = ~isnan(gamma_tr);

if any(i1)
    semilogx(alpha_grid(i1),gamma_fdi(i1),'o-','LineWidth',1.0, ...
        'DisplayName','Full-direction invariance');
end

if any(i2)
    semilogx(alpha_grid(i2),gamma_tr(i2),'s-','LineWidth',1.0, ...
        'DisplayName','Structure transport');
end

grid on;
xlabel('\alpha');
ylabel('\gamma');
legend('show','Location','best');
function out = solve_slice(a,isTR,A,B2,B1,C,D,opt)
yalmip('clear');
P=sdpvar(6,6,'symmetric'); Z=sdpvar(6,6,'symmetric');
X0=sdpvar(3,3,'full'); R0=sdpvar(2,3,'full'); g=sdpvar(1);
if isTR
    X=[X0 -X0; zeros(3) X0];
    R=[R0 zeros(2,3); zeros(2,3) R0];
else
    X=[X0 zeros(3); zeros(3) X0];
    R=[R0 R0; zeros(2,3) R0];
end
Q=[A*X-B2*R; -X; C*X-D*R]*[eye(6) a*eye(6) zeros(6,10)];
M1=-[zeros(6) P zeros(6,10); P zeros(6) zeros(6,10); ...
     zeros(10,6) zeros(10,6) -eye(10)]-Q-Q';
M2=[Z B1'; B1 P];
sol=optimize([M1>=0,M2>=0,P>=0,Z>=0,trace(Z)<=g],g,opt);
out=struct('ok',false,'gamma',inf,'h2sq',nan,'X',[],'R',[],'K',[]);
if sol.problem~=0, return; end
Xv=value(X); Rv=value(R); gv=value(g);
if ~isfinite(gv) || rcond(Xv)<1e-10, return; end
Kv=Rv/Xv;
Acl=A-B2*Kv;
if max(real(eig(Acl)))>=0, return; end
% Actual closed-loop H2^2 from w to z.
% Closed loop: xdot = Acl*x + B1*w, z = (C-D*Kv)*x.
Ccl=C-D*Kv;
sys = ss(Acl,B1,Ccl,zeros(10,6));
h2sq=norm(sys,2)^2;
if ~isfinite(h2sq) || h2sq < -1e-8, return; end
h2sq=max(h2sq,0);
out=struct('ok',true,'gamma',gv,'h2sq',h2sq,'X',Xv,'R',Rv,'K',Kv);
end


%%%%%% realized h2 cost 


K_co_ref=[1.7962, 2.7236, -1.5147; 
    0.1511, 0.1290, 5.6098];
K_ref=[K_co_ref K_co_ref;zeros(2,3) K_co_ref];
Acl_ref=A-B2*K_ref;
Ccl_ref=C-D*K_ref;
sys_ref = ss(Acl_ref,B1,Ccl_ref,zeros(10,6));
h2sq_ref = norm(sys_ref,2)^2

K_co_pro=[3.1683,   4.2569,   -1.4781;
    1.0851,    0.7454,    2.9678];
K_pro=[K_co_pro K_co_pro;zeros(2,3) K_co_pro];

Acl_pro=A-B2*K_pro;
Ccl_pro=C-D*K_pro;
sys_pro = ss(Acl_pro,B1,Ccl_pro,zeros(10,6));
h2sq_ref = norm(sys_pro,2)^2

eig(Acl_pro)