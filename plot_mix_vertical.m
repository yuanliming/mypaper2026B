clear; close all; clc;

files = {'tau0.2.mat','tau0.5.mat'};
taus  = [0.2 0.5];
alpha = logspace(-2,2,101);

labels = {'$F_1$','$F_2$','$F_{12}^{+}$','$F_{12}^{-}$', ...
          '$F_{23}^{+}$','$F_{23}^{-}$'};

mk = {'o','s','d','^','v','>'};
ls = {'-','--','-.',':','-','--'};

figure('Units','centimeters','Position',[2 2 9 13]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

for k = 1:2
    S = load(files{k});
    G = S.gamma_all;

    nexttile; hold on;

    for i = 1:6
        g = G(i,:);

        % Keep main continuous feasible branch
        j = find(~isfinite(g),1);
        if isempty(j), j = numel(g)+1; end
        idx = 1:j-1;

        h(i) = semilogx(alpha(idx),g(idx),[ls{i} mk{i}], ...
            'LineWidth',1.0, ...
            'MarkerSize',4, ...
            'MarkerIndices',1:6:numel(idx));
    end

    grid on;
    box on;

    % Emphasize the neighborhood of the best bounds
    xlim([2e-2 2]);
    ylim([2 3.2]);

    title(sprintf('\\tau = %.1f',taus(k)));
    ylabel('$\gamma$','Interpreter','latex');
    if k == 2
        xlabel('$\alpha$','Interpreter','latex');
    end
end

lgd = legend(h,labels, ...
    'Interpreter','latex', ...
    'Orientation','horizontal', ...
    'NumColumns',3);

lgd.Layout.Tile = 'south';

exportgraphics(gcf,'gamma_alpha_tau_zoom.pdf', ...
    'ContentType','vector');

exportgraphics(gcf,'gamma_alpha_tau_zoom.png', ...
    'Resolution',300);