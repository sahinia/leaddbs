function fnames = ea_sweetspot_export_empirical_sigmap(permtestFile, alphalevel, posvisible, negvisible, outdir)
% Exports the empirical R-map from a saved permtest .mat file (group 1
% only, matching predpermtest's own convention), masked to voxels
% significant at pemp <= alphalevel -- this reproduces exactly the
% significance mask ea_sweetspot's predpermtest applies to the empirical
% training map before out-of-sample prediction (see its Remp_train/nonsig
% logic), plus the positive/negative visibility mask applied at prediction
% time (see ea_sweetspot's private maskvals: posvisible=0 blanks R>0,
% negvisible=0 blanks R<0).
%
% Works against either a <ID>_permtest.mat file or its _slim.mat companion
% -- Remp/pemp are carried over unchanged in the slim file, so both have
% what this needs.
%
% permtestFile - path to the training permtest (or slim) .mat file.
% alphalevel   - uncorrected significance cutoff on pemp (default 0.05).
% posvisible   - include positive ("sweet spot") voxels (default 1).
% negvisible   - include negative ("sour spot") voxels (default 1).
% outdir       - output folder (default: same folder as permtestFile).
%
% Writes <basename>_empirical_sig_p<alpha>_r/_l.nii and returns the
% written file paths.

if ~exist('alphalevel', 'var') || isempty(alphalevel)
    alphalevel = 0.05;
end
if ~exist('posvisible', 'var') || isempty(posvisible)
    posvisible = 1;
end
if ~exist('negvisible', 'var') || isempty(negvisible)
    negvisible = 1;
end
if ~exist('outdir', 'var') || isempty(outdir)
    outdir = [fileparts(permtestFile), filesep];
end

S = load(permtestFile, 'Remp', 'pemp', 'space');
nsides = numel(S.space);
[~, basefname] = fileparts(permtestFile);

sigMap = cell(1, nsides);
for side = 1:nsides
    Remp = S.Remp{1, side};
    pemp = S.pemp{1, side};

    nonsig = isnan(pemp) | pemp > alphalevel;
    v = Remp;
    v(nonsig) = nan;

    if ~posvisible
        v(v > 0) = nan;
    end
    if ~negvisible
        v(v < 0) = nan;
    end

    sigMap{side} = v;
    fprintf('Side %d: %d/%d voxels significant at p<=%.3g (posvisible=%d, negvisible=%d).\n', ...
        side, sum(~isnan(v)), numel(v), alphalevel, posvisible, negvisible);
end

fnames = ea_sweetspot_vals2nii(S.space, sigMap, outdir, sprintf('%s_empirical_sig_p%.3g', basefname, alphalevel));
end
