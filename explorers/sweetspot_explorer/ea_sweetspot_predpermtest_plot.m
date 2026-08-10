function ea_sweetspot_predpermtest_plot(predpermtestFile, responselabel)
% Regenerates the two prediction-validation figures (empirical correlation
% scatter + out-of-sample null distribution) from a saved
% <ID>_predpermtest.mat file, without rerunning the prediction. Standalone --
% only needs the file path, no live ea_sweetspot object.
%
% responselabel (optional): x-axis label for the correlation scatter plot.
% Defaults to the label saved at prediction time, or 'Response Variable' if
% that's unavailable (files saved before that field was added).

pred = load(predpermtestFile);

if ~exist('responselabel', 'var') || isempty(responselabel)
    if isfield(pred, 'responsevarlabel') && ~isempty(pred.responsevarlabel)
        responselabel = pred.responsevarlabel;
    else
        responselabel = 'Response Variable';
    end
end

% Figure 1: empirical Ihat vs. the test cohort's real, unpermuted responsevar.
% Iemp/Ihatemp were only added to predpermtest's saved output recently -- older
% files won't have them. Skip this figure rather than erroring the whole
% function, since Figure 2 below doesn't depend on either field.
if isfield(pred, 'Iemp') && isfield(pred, 'Ihatemp')
    h1 = ea_corrbox(pred.Iemp, pred.Ihatemp, pred.pperm, ...
        {sprintf('Out-of-Sample Prediction (%s)', pred.testID), responselabel, 'Predicted Score (Ihat)'});
    try % ea_corrbox is a shared, gramm-based utility -- best-effort cleanup, non-fatal if its internal structure doesn't cooperate
        set(h1, 'Color', 'w');
        axesInH1 = findall(h1, 'Type', 'axes');
        for a = 1:numel(axesInH1)
            set(axesInH1(a), 'Color', 'w');
            box(axesInH1(a), 'off');
        end
    end
else
    warning('This predpermtest.mat file was saved before Iemp/Ihatemp were added -- skipping the correlation scatter (Figure 1). Rerun predpermtest to regenerate a file with this figure available. Showing the null-distribution plot (Figure 2) only.');
end

% Figure 2: null distribution of out-of-sample R, empirical R marked, rank
% annotated. tail must match whatever predpermtest actually used to compute
% pred.pperm/exceedCount (Figure 1's ea_corrbox call above uses pred.pperm
% directly, so it's already tail-correct; this plot instead recomputes its
% own rank/p-value from Rpredperm/Rpredemp, so passing the wrong tail here
% would show a DIFFERENT p-value than Figure 1 for the same file).
% pred.tail only exists on files saved after the tail option was added --
% older files were always two-sided, so 'both' is the correct fallback for
% them, not a guess.
if isfield(pred, 'tail') && ~isempty(pred.tail)
    plotTail = pred.tail;
else
    plotTail = 'both';
end
ea_sweetspot_nulldist_plot(pred.Rpredperm, pred.Rpredemp, ...
    sprintf('Out-of-sample prediction R (%s)', pred.corrtype), ...
    sprintf('Out-of-Sample Permutation Test (%d/%d permutations produced no defined prediction)', pred.nNaNperm, numel(pred.Rpredperm)), ...
    plotTail);
