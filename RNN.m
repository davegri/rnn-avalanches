tic
clc; clear;
sigmas = [0.7, 0.8, 0.9, 0.95, 1, 1.05, 1.1];
for s = 1:length(sigmas)
    aval{s} = binary_RNN(sigmas(s));
end
toc

%% kappa
for s = 1:length(sigmas)
    kappas(s) = get_kappa(aval{s});
end
plot(sigmas, kappas, '-*', 'LineWidth', 2);
title("\kappa(\sigma) in binary RNN simulation", "FontSize", 16);
xlabel("\sigma", "FontSize", 18, 'FontWeight','bold');
ylabel("\kappa", "FontSize", 18, 'FontWeight','bold');

hold on;
plot(sigmas,sigmas, 'LineWidth', 2);
legend("\kappa(\sigma)", "\kappa = \sigma", 'Location', 'NorthWest', 'AutoUpdate','off');
yline(1)
xline(1)

%% plot
tiledlayout('flow');
colormap summer
% pdf
s = 1:100:10e8;
Ps = s.^(-3/2);
Ps = Ps./sum(Ps);
nexttile;
c_map = winter(length(sigmas));
for i = 1:length(sigmas)
    color = c_map(i,:);
    if sigmas(i) == 1.0
        color = 'red';
    end
    plot_avalanches_pdf(aval{i}, color);
    hold on;
end
title("Avalanche size pdf for various branching parameters", "FontSize", 16);
plot(s,Ps, '--', 'LineWidth', 1);
leg = legend(["\sigma = " + string(sigmas) "s^{-3/2}"], 'FontSize', 14, 'Location', 'SouthWest');
xlabel("Avalanche Size", 'FontSize', 14);
ylabel("Probability Density", 'FontSize', 14);

%%
%cdf
s = 1:1:10e5;
Ps = s.^(-3/2);
Ps = Ps./sum(Ps);
nexttile([1 3]);
for j = 1:length(sigmas)
    color = c_map(j,:);
    if sigmas(j) == 1.0
        color = 'red';
    end
    plot_avalanches_cdf(aval{j}, color);
    hold on;
end

s_cdf = cumsum(Ps);
plot(s,s_cdf, '--', 'LineWidth', 2);
sgtitle("CDF's for various branching parameters", "FontSize", 16);
ylim([min(s_cdf) 1.02]);
legend(["\sigma = " + string(sigmas) "s^{-3/2}"], 'Location', 'SouthEast', 'FontSize', 14);
xlabel("Avalanche Size", 'FontSize', 14);
ylabel("Probability", 'FontSize', 14);

function plot_avalanches_pdf(avalanches, color)
    edges = logspace(0,9,20);
    N = histcounts(avalanches, edges, 'Normalization','pdf');
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    loglog(edges, N, 'LineWidth', 1.5, 'Color', color);
end
function plot_avalanches_cdf(avalanches, color)
    edges = logspace(0, 9,30);
    N = histcounts(avalanches, edges, 'Normalization','cdf');
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges, N, 'LineWidth', 1, 'Color', color);
    ax = gca;
    set(ax, "XScale", "log");
    set(ax, "YScale", "linear")
end

function kappa = get_kappa(aval)
    m = 10;
    edges = logspace(0, 3 ,20);
    N = histcounts(aval, edges, 'Normalization','cdf');
    edges = edges(2:end);
    l = min(edges);
    L = max(edges);
    betas = logspace(l, log10(L), m);

    for i=1:length(betas)
        F_data(i) = N(find(edges<=betas(i),1,'last'));
    end

    F_theo = (1-sqrt(l./betas))/(1-sqrt(l/L));
    kappa = 1 + mean(F_theo-F_data);
end

function avalanche_size = binary_RNN(sigma)
    num_runs = 10000;
    N = 500;
    avalanche_size = zeros(num_runs,1);
    p = 1/N*rand(N);
    p = p./sum(p,2); % normalize each row of p so that it sums to 1
    p = sigma*p; % multiply by branching factor

    parfor r = 1:num_runs
        s = zeros(N,1);
        s(randi(N)) = 1;
        aval = 1;
        while any(s)
            zeta = rand(N,1);
            p_J = 1 - prod((1-p(:,s==1)),2);
            s = p_J > zeta;
            aval = aval + sum(s);
%             if aval > 500000
%                 break
%             end
        end
        avalanche_size(r) = aval;
    end
end
