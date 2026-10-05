%% example2_mixed_affine.m
% Example 2: mixed affine-equality H2 synthesis.
% Reproduces the candidate screening, Table results, and Fig. 1.
%
% Controller constraints:
%   k11 = 0,  k23 = 1,  k12 = k21,  k22 + 2*k13 = 3.
%
% Requires YALMIP, SDPT3, and Control System Toolbox.

clear; clc; close all;

alpha_grid = logspace(-2,2,101);
opt = sdpsettings('solver','sdpt3','verbose',0);

%% Plant and affine controller set
A  = [-1 2 1; 2 -1 0; 4 4 -2];
B2 = [2 0; 0 2; 2 2];
B1 = eye(3);
C  = [eye(3); zeros(2,3)];
D  = [zeros(3,2); eye(2)];

Kp = [0 0 0; 0 3 1];
N1 = [0 1 0; 1 0 0];
N2 = [0 0 1; 0 -2 0];

fprintf('Open-loop spectral abscissa = %.4f\n',max(real(eig(A))));

%% Candidate library in (56)
Fb = diag([1 1 2]);
E12 = [0 1 0; 0 0 0; 0 0 0];
E23 = [0 0 0; 0 0 1; 0 0 0];

F = cell(10,1);
name = cell(10,1);
tau = nan(10,1);

F{1}=eye(3); F{2}=Fb;
name{1}='F1 = I3'; name{2}='F2 = Fb';

k = 2;
for t = [0.2 0.5]
    F{k+1}=Fb*(eye(3)+t*E12); name{k+1}='F12+'; tau(k+1)=t;
    F{k+2}=Fb*(eye(3)-t*E12); name{k+2}='F12-'; tau(k+2)=t;
    F{k+3}=Fb*(eye(3)+t*E23); name{k+3}='F23+'; tau(k+3)=t;
    F{k+4}=Fb*(eye(3)-t*E23); name{k+4}='F23-'; tau(k+4)=t;
    k = k+4;
end

%% Screening and H2 synthesis
nF = numel(F);
gamma_all = nan(nF,numel(alpha_grid));
h2_all = nan(nF,numel(alpha_grid));
best = repmat(struct('ok',false,'gamma',inf,'h2',nan,'alpha',nan,'K',[]),nF,1);

fprintf('\nScreening:\n');
for j = 1:nF
    % Here A_H = span{I3}; He(F)>0 directly gives denominator compatibility.
    den_ok = min(eig(F{j}+F{j}')) > 1e-10;
    cert_status = certificate_test(F{j},A,B2,Kp,N1,N2,opt);

    if cert_status == 1
        cert_msg = 'INFEASIBLE -> survives';
    elseif cert_status == 0
        cert_msg = 'FEASIBLE -> discarded';
    else
        cert_msg = yalmiperror(cert_status);
    end

    fprintf('  %-8s',name{j});
    if ~isnan(tau(j)), fprintf(' tau=%.1f',tau(j)); else, fprintf('        '); end
    fprintf(': denominator %s, certificate %s\n',passfail(den_ok),cert_msg);

    if ~den_ok || cert_status ~= 1
        continue
    end

    for ia = 1:numel(alpha_grid)
        s = solve_slice(alpha_grid(ia),F{j},A,B2,B1,C,D,Kp,N1,N2,opt);
        if s.ok
            gamma_all(j,ia) = s.gamma;
            h2_all(j,ia) = s.h2;
            if s.gamma < best(j).gamma
                best(j) = s;
                best(j).alpha = alpha_grid(ia);
            end
        end
    end
end

%% Table values: best tested member of each reference family
groups = {1, 2, [3 7], [4 8], [5 9], [6 10]};
gname  = {'F1 = I3 (self)','F2 = Fb','F12+','F12-','F23+','F23-'};

fprintf('\nMixed-affine H2 synthesis results\n');
fprintf('%-16s %6s %10s %10s %10s\n','Reference','tau','alpha','gamma','H2^2');

family_best = nan(6,1);
for g = 1:6
    ids = groups{g};
    [~,q] = min(arrayfun(@(i)best(i).gamma,ids));
    j = ids(q);
    family_best(g) = j;

    if isnan(tau(j)), ts='--'; else, ts=sprintf('%.1f',tau(j)); end
    fprintf('%-16s %6s %10.4f %10.4f %10.4f\n', ...
        gname{g},ts,best(j).alpha,best(j).gamma,best(j).h2);
end

%% Best certified bound and smallest realized cost among family winners
[~,g_cert] = min(arrayfun(@(j)best(j).gamma,family_best));
j_cert = family_best(g_cert);

[~,g_real] = min(arrayfun(@(j)best(j).h2,family_best));
j_real = family_best(g_real);

fprintf('\nBest certified bound:\n');
fprintf('  %s, tau=%.1f, alpha=%.4f, gamma=%.4f, H2^2=%.4f\n', ...
    name{j_cert},tau(j_cert),best(j_cert).alpha,best(j_cert).gamma,best(j_cert).h2);
disp('  K ='); disp(best(j_cert).K);

improvement = 100*(best(1).gamma-best(j_cert).gamma)/best(1).gamma;
fprintf('  bound reduction relative to self transport = %.1f%%\n',improvement);

fprintf('\nSmallest realized cost among family winners:\n');
fprintf('  %s, tau=%.1f, H2^2=%.4f\n', ...
    name{j_real},tau(j_real),best(j_real).h2);

%% Fig. 1: certified bounds for tau = 0.2 and 0.5
% Use color + different line/marker combinations so the figure remains
% distinguishable in grayscale printing.
plot_ids = {[1 2 3 4 5 6],[1 2 7 8 9 10]};
plot_labels = {'$F_1$','$F_2$','$F_{12}^{+}$','$F_{12}^{-}$', ...
               '$F_{23}^{+}$','$F_{23}^{-}$'};
tau_plot = [0.2 0.5];

line_style = {'-','--','-.',':','-','--'};
marker = {'o','s','d','^','v','>'};

figure('Units','centimeters','Position',[2 2 9 13],'Color','w');
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
h = gobjects(6,1);

for p = 1:2
    ax = nexttile; hold(ax,'on');
    ids = plot_ids{p};

    for r = 1:6
        j = ids(r);
        ok = isfinite(gamma_all(j,:));
        npt = nnz(ok);
        h(r) = semilogx(alpha_grid(ok),gamma_all(j,ok), ...
            'LineStyle',line_style{r}, ...
            'Marker',marker{r}, ...
            'MarkerIndices',1:6:npt, ...
            'LineWidth',1.0,'MarkerSize',4);
    end

    grid on; box on;
    xlim([2e-2 2]); ylim([2 3.2]);
    set(ax,'FontName','Times New Roman','FontSize',8,'LineWidth',0.8);

    title(sprintf('$\\tau = %.1f$',tau_plot(p)),'Interpreter','latex');
    ylabel('$\gamma$','Interpreter','latex');
    if p==2
        xlabel('$\alpha$','Interpreter','latex');
    end
end

lgd = legend(h,plot_labels, ...
    'Interpreter','latex', ...
    'Orientation','horizontal', ...
    'NumColumns',3, ...
    'FontSize',8);
lgd.Layout.Tile = 'south';

% % Vector PDF for the manuscript and 300-dpi PNG for quick viewing.
% exportgraphics(gcf,'example2_gamma_alpha.pdf','ContentType','vector');
% exportgraphics(gcf,'example2_gamma_alpha.png','Resolution',300);

%% ------------------------------------------------------------------------
function status = certificate_test(F,A,B2,Kp,N1,N2,opt)
% Recovery-core certificate SDP for A_H = span{I3}, S_H = span{N1,N2}.
yalmip('clear');
Q = sdpvar(3,3,'symmetric');
Y = sdpvar(3,3,'symmetric');
Ap = A-B2*Kp;

BF = B2'*Q*F';
AF = (Ap'*Q-Y)*F';

con = [Q>=0, Y>=0, trace(Q)==1, ...
       sum(sum(BF.*N1))==0, sum(sum(BF.*N2))==0, trace(AF)==0];

sol = optimize(con,[],opt);
status = sol.problem;
end

function out = solve_slice(a,F,A,B2,B1,C,D,Kp,N1,N2,opt)
% Fixed-slice parameterization:
% X = xi F,  R = (xi Kp + c1 N1 + c2 N2)F.
yalmip('clear');

P = sdpvar(3,3,'symmetric');
Z = sdpvar(3,3,'symmetric');
xi = sdpvar(1);
c = sdpvar(2,1);
g = sdpvar(1);

X = xi*F;
R = (xi*Kp + c(1)*N1 + c(2)*N2)*F;

Qm = [A*X-B2*R; -X; C*X-D*R] * ...
     [eye(3) a*eye(3) zeros(3,5)];

M0 = [zeros(3) P zeros(3,5);
      P zeros(3) zeros(3,5);
      zeros(5,3) zeros(5,3) -eye(5)];

L1 = -M0-Qm-Qm';
L2 = [Z B1'; B1 P];

sol = optimize([L1>=1e-8*eye(11), L2>=1e-8*eye(6), ...
                P>=1e-8*eye(3), Z>=1e-8*eye(3), trace(Z)<=g],g,opt);

out = struct('ok',false,'gamma',inf,'h2',nan,'alpha',nan,'K',[]);
if sol.problem~=0, return; end

Xv=value(X); Rv=value(R);
if rcond(Xv)<1e-10, return; end

K=Rv/Xv;
Acl=A-B2*K;
if max(real(eig(Acl)))>=-1e-8, return; end

sys=ss(Acl,B1,C-D*K,zeros(5,3));
out=struct('ok',true,'gamma',value(g),'h2',norm(sys,2)^2, ...
           'alpha',nan,'K',K);
end

function s = passfail(tf)
if tf, s='PASS'; else, s='FAIL'; end
end
