function fig = benchmarkComparisonPlot(results)
% Plot comparisons between libtopotoolbox and native MATLAB implementations
%
% 1. Run benchmarks with `buildtool benchmark`
% 2. Load results from the generated file with `results =
% readtable("tests/results/benchmark.csv","TextType","string")`
% 3. Plot with `benchmarkComparisonPlot(results)`
%
% The resulting figure has one plot for each snapshot test method.
% Each plot shows the median benchmark result for the snapshot test using
% libtopotoolbox (y-axis) against the snapshot test using the native MATLAB
% implementation (x-axis). Each snapshot dataset is plotted as a line with
% '+' markers, and the most recent benchmark run is plotted as a filled,
% colored circle. A zone of +-10% of the baseline time is shaded gray.

% Filter out groups that don't have libtopotoolbox benchmarks
G = groupfilter(results,["Site","Function","RunIdentifier"],@(x) nnz(x) > 0,"UseLibTT");
sites = unique(G.Site);
functions = unique(G.Function);


% Comparison plot

% Extract run information
Runs = pivot(G,Rows="RunIdentifier",DataVariable="RunDate",Method="min");

% Pivot results table to expose the median timings
P = pivot(G,Columns=["Function","Site","UseLibTT"],Rows="RunIdentifier",DataVariable="Median");

% Sort pivot table by the run date
P = join(P,Runs);
P = sortrows(P,"min_RunDate");

fig = figure('WindowStyle','docked');
t = tiledlayout("flow");
for f = functions'
    nexttile;
    for s = sites'
        r = P.(f).(s);
        loglog(r,"0","1",'marker','+','Color','k','HandleVisibility','off');
        hold on
        scatter(r(end,:),"0","1","filled",'DisplayName',s);
    end

    % Plot the shaded zone from y=0.9x to y = 1.1x
    xs = xlim;
    ys = ylim;
    xmin = 0.9*min(xs(1),ys(1));
    xmax = 1.1*max(xs(2),ys(2));
    fill([xmin;xmin;xmax;xmax], ...
        [0.9*xmin;1.1*xmin;1.1*xmax;0.9*xmax], ...
        [0.5,0.5,0.5], ...
        'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off');
    hold off
    xlabel("Baseline (s)");
    ylabel('');
    title(f);
    xlim([xmin,xmax]);
    ylim([xmin,xmax]);
    axis square manual;
end
ylabel(t,"Libtopotoolbox (s)");
lgd = legend('Interpreter','none');
lgd.Layout.Tile = 'east';

end