function [h, pval, exceedCount] = ea_sweetspot_nulldist_plot(nullvals, empval, xlabelstr, titlestr, tail)
% Plots a permutation null distribution with the empirical value marked
% and its rank/p-value annotated in the title. Publication style: white
% background, no box, thin lines. Shared by the voxelwise-map-level plots
% (omnibus, max-statistic) and the out-of-sample prediction plot.
%
% tail: 'right' (nullvals >= empval, e.g. omnibus Q, positive max-stat)
%       'left'  (nullvals <= empval, e.g. negative max-stat)
%       'both'  (abs(nullvals) >= abs(empval), e.g. out-of-sample prediction R)
% Defaults to 'right'.
%
% NaN entries in nullvals are never counted as exceeding (comparisons
% against NaN are always false), but remain in the denominator -- a
% permutation that produced no defined statistic is evidence against the
% null, not missing data (see ea_sweetspot.m predpermtest for the same
% convention/reasoning).

if ~exist('tail', 'var') || isempty(tail)
    tail = 'right';
end

pal = ea_sweetspot_palette();

n = numel(nullvals);
switch tail
    case 'right'
        exceedCount = sum(nullvals >= empval);
    case 'left'
        exceedCount = sum(nullvals <= empval);
    case 'both'
        exceedCount = sum(abs(nullvals) >= abs(empval));
end
pval = exceedCount / n;

h = figure('Color', 'w');
histogram(nullvals, 'FaceColor', pal.blue, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
hold on;
xline(empval, 'Color', pal.red, 'LineWidth', 2);
legend({'Null distribution', 'Empirical'}, 'Box', 'off', 'Location', 'best');
xlabel(xlabelstr);
ylabel('Number of permutations');
title({titlestr, sprintf('Empirical = %.3f, ranked %d of %d (p = %.3f)', empval, exceedCount, n, pval)});
ea_sweetspot_style_axes(gca);
