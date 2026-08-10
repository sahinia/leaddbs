function results = ea_sweetspot_permtest_stats(permtestFile, alphaVoxelwise, alphaEisenstein, alphaMaxstat, useTailApprox, batchSize)
% Reads a <ID>_permtest.mat file (produced by ea_sweetspot's permtest method)
% and computes three complementary permutation-based significance analyses,
% entirely standalone -- only the file path and explicit alpha levels are
% needed, no live ea_sweetspot object or GUI settings.
%
% alphaVoxelwise  - per-voxel uncorrected permutation threshold (each voxel's
%                    empirical value ranked against its OWN permuted null,
%                    independently per voxel; no cross-voxel correction).
%                    Default 0.05.
% alphaEisenstein - per-voxel significance cutoff used to decide which voxels
%                    are summed into the Eisenstein et al. 2014 (Annals of
%                    Neurology 76:279-295) omnibus Q statistic:
%                    Q = sum(-log10(p)) over voxels with p <= alphaEisenstein,
%                    computed for the empirical map and every permuted map,
%                    giving a permutation-based whole-map significance test.
%                    Default 0.05.
% alphaMaxstat    - tail probability for the max-statistic FWER correction
%                    (e.g. 0.05 -> the 95th percentile of the per-permutation
%                    maximum, 5th percentile of the per-permutation minimum).
%                    Default 0.05.
% useTailApprox   - if true, the Eisenstein and max-statistic p-values (each
%                    a single test against one shared null, unlike the
%                    per-voxel test) are refined with a Generalized Pareto
%                    Distribution tail fit when the rank-based p is already
%                    small (see ea_sweetspot_gpd_pval, ported from FSL
%                    PALM's palm_pareto.m) -- lets you report p below the
%                    1/(Nperm+1) permutation-count floor when the null's
%                    tail actually fits a GPD well. Default false (matches
%                    prior behavior exactly).
% batchSize       - for group 1 (see below), Rperm/pperm are read off disk
%                    in row batches of this many permutations at a time
%                    rather than loaded whole. Default 250; lower it
%                    further on tighter-memory machines, at the cost of
%                    more, smaller disk reads.
%
% Positive ("sweet spot") and negative ("sour spot") directions are tracked
% and thresholded separately throughout, rather than folded into |R|.
%
% Voxelwise-threshold and max-statistic maps are exported as nifti files next
% to permtestFile (one file per side, '_r'/'_l' suffixed, matching how
% permtest itself exports Remp). The Eisenstein and max-statistic null
% distributions are plotted with the empirical value marked.
%
% Returns a struct with everything computed (per group/side, matching the
% {group,side} cell convention used throughout the sweetspot explorer).

if ~exist('alphaVoxelwise', 'var') || isempty(alphaVoxelwise)
    alphaVoxelwise = 0.05;
end
if ~exist('alphaEisenstein', 'var') || isempty(alphaEisenstein)
    alphaEisenstein = 0.05;
end
if ~exist('alphaMaxstat', 'var') || isempty(alphaMaxstat)
    alphaMaxstat = 0.05;
end
if ~exist('useTailApprox', 'var') || isempty(useTailApprox)
    useTailApprox = false;
end
if ~exist('batchSize', 'var') || isempty(batchSize)
    batchSize = 250;
end

% Dot-free string forms of the alpha values for use in exported FILENAMES
% only (results.alphaVoxelwise/alphaMaxstat above keep the real numeric
% value) -- ea_write_nii's ea_stripext call treats the FIRST '.' anywhere
% in a filename as an extension boundary, not just the one before '.nii',
% so a name like '..._maxstat_p0.05_group1.nii' silently gets written to
% disk as '..._maxstat_p0.nii' (losing the '05_group1' part entirely, and
% colliding with any other export whose name also truncates to the same
% thing -- e.g. a different group, or a different alpha). Stripping the
% dot instead of embedding it sidesteps that entirely.
alphaVoxelwiseStr = strrep(sprintf('%.3g', alphaVoxelwise), '.', '');
alphaMaxstatStr = strrep(sprintf('%.3g', alphaMaxstat), '.', '');

% Only the small variables (Nvoxels x 1 or smaller) are loaded whole here.
% Rperm/pperm (Nperm x Nvoxels per side) are the memory-heavy part of a
% saved permtest file and are read separately, in row-chunked batches, via
% ea_sweetspot_permtest_slim + matfile() below -- a blanket load() of
% those used to push MATLAB to ~159GB and effectively hang the machine on
% large (high-Nperm, high-resolution) runs.
train = load(permtestFile, 'space', 'Remp', 'pemp', 'Nperm');
outdir = [fileparts(permtestFile), filesep];
[~, basefname] = fileparts(permtestFile);

results = struct();
results.sourcefile = permtestFile;
results.alphaVoxelwise = alphaVoxelwise;
results.alphaEisenstein = alphaEisenstein;
results.alphaMaxstat = alphaMaxstat;
results.useTailApprox = useTailApprox;

nsides = numel(train.space);
Nperm = train.Nperm;

% Group 1 is always readable this way: permtest() now saves it directly in
% this flat, single-precision, chunk-readable layout; ea_sweetspot_permtest_slim
% converts an older/legacy file to match, or passes an already-flat one
% through unchanged (see that function). Groups beyond 1 (obj.splitbygroup)
% have never been supported by this flat layout -- trainHeavy, loaded lazily
% only if a group > 1 is actually encountered below, keeps that rare path
% working exactly as before.
slimFile = ea_sweetspot_permtest_slim(permtestFile);
mTrain = matfile(slimFile);
trainHeavy = [];

summaryLines = {}; % one line per group/side, appended to this run's README at the end

for gi = 1:size(train.Remp, 1)
    if all(cellfun(@isempty, train.Remp(gi,:)))
        continue
    end

    voxThreshMap = cell(1, nsides);
    maxstatMap = cell(1, nsides);

    for side = 1:nsides
        Remp = train.Remp{gi,side};   % Nvoxels x 1
        pemp = train.pemp{gi,side};   % Nvoxels x 1

        tag = sprintf('group%d side%d', gi, side);

        %% 1+3 setup: accumulate the per-voxel exceedance counts (#1) and the
        % per-permutation Eisenstein Q / max / min (#2, #3) in row batches,
        % rather than ever holding the full Nperm x Nvoxels Rperm/pperm in
        % memory at once.
        RempRow = Remp'; % 1 x Nvoxels, for broadcasting against Rperm's rows
        countPos = zeros(1, numel(Remp));
        countNeg = zeros(1, numel(Remp));
        Qperm = nan(Nperm, 1);
        maxRperm = nan(Nperm, 1);
        minRperm = nan(Nperm, 1);

        if gi == 1
            for batchStart = 1:batchSize:Nperm
                batchIdx = batchStart:min(batchStart+batchSize-1, Nperm);
                Rbatch = double(mTrain.(sprintf('Rperm_side%d', side))(batchIdx,:));
                pbatch = double(mTrain.(sprintf('pperm_side%d', side))(batchIdx,:));

                countPos = countPos + sum(Rbatch >= RempRow, 1);
                countNeg = countNeg + sum(Rbatch <= RempRow, 1);

                sigMaskBatch = pbatch <= alphaEisenstein;
                logpBatch = -log10(max(pbatch, eps));
                logpBatch(~sigMaskBatch) = 0;
                Qperm(batchIdx) = ea_nansum(logpBatch, 2);

                % NOTE: ea_nanmax/ea_nanmin (ext_libs/nan) do NOT use MATLAB's
                % own max(A,[],dim) convention -- their 3-argument form is
                % (a,dim,b) for an ELEMENTWISE max/min of two same-sized
                % arrays when dim is empty, not "reduce along dim". The
                % 2-argument form (a,dim) is what reduces along a dimension.
                maxRperm(batchIdx) = ea_nanmax(Rbatch, 2);
                minRperm(batchIdx) = ea_nanmin(Rbatch, 2);
            end
        else
            % Rare multi-group path -- Rperm/pperm only exist in the
            % original nested-cell, double-precision form for group > 1.
            if isempty(trainHeavy)
                trainHeavy = load(permtestFile, 'Rperm', 'pperm');
            end
            Rperm = trainHeavy.Rperm{gi,side};
            pperm = trainHeavy.pperm{gi,side};

            countPos = sum(Rperm >= RempRow, 1);
            countNeg = sum(Rperm <= RempRow, 1);

            sigMaskPerm = pperm <= alphaEisenstein;
            logpPerm = -log10(max(pperm, eps));
            logpPerm(~sigMaskPerm) = 0;
            Qperm = ea_nansum(logpPerm, 2);

            maxRperm = ea_nanmax(Rperm, 2);
            minRperm = ea_nanmin(Rperm, 2);
        end

        % +1 in numerator and denominator: the empirical map is itself one
        % valid draw under the null, so it belongs in its own reference set.
        % Without it, a voxel with zero exceedances would report p = 0,
        % which is not a valid claim from a finite number of permutations
        % (Phipson & Smyth 2010).
        pVoxPos = (countPos + 1) / (Nperm + 1);
        pVoxNeg = (countNeg + 1) / (Nperm + 1);
        pVoxPos(isnan(RempRow)) = NaN;
        pVoxNeg(isnan(RempRow)) = NaN;

        voxSurvivePos = pVoxPos <= alphaVoxelwise;
        voxSurviveNeg = pVoxNeg <= alphaVoxelwise;

        thisVoxMap = nan(size(Remp));
        thisVoxMap(voxSurvivePos') = Remp(voxSurvivePos');
        thisVoxMap(voxSurviveNeg') = Remp(voxSurviveNeg');
        voxThreshMap{side} = thisVoxMap;

        results.voxelwise.pPos{gi,side} = pVoxPos;
        results.voxelwise.pNeg{gi,side} = pVoxNeg;
        results.voxelwise.thresholdedMap{gi,side} = thisVoxMap;

        %% 2. Eisenstein 2014 omnibus Q statistic: sum(-log10 p) over per-voxel-significant voxels
        sigMaskEmp = pemp' <= alphaEisenstein;
        logpEmp = -log10(max(pemp', eps));
        logpEmp(~sigMaskEmp) = 0;
        Qemp = ea_nansum(logpEmp);

        results.eisenstein.Qemp{gi,side} = Qemp;
        results.eisenstein.Qperm{gi,side} = Qperm;

        [~, results.eisenstein.p{gi,side}, results.eisenstein.rank{gi,side}] = ea_sweetspot_nulldist_plot( ...
            Qperm, Qemp, 'Q = \Sigma -log_{10}(p)', ...
            sprintf('Eisenstein 2014 Omnibus Test (%s, \\alpha=%.3g)', tag, alphaEisenstein), 'right', useTailApprox);

        %% 3. Max-statistic FWER correction (pos/neg separate)
        maxRemp = ea_nanmax(Remp);
        minRemp = ea_nanmin(Remp);

        % Exact-rank threshold, consistent with the (count+1)/(Nperm+1)
        % p-value formula used above and in ea_sweetspot_nulldist_plot --
        % NOT prctile(), whose interpolated quantile can land on a value
        % that never actually occurred in Rperm, so the exported thresholded
        % map and the printed p-value for the same test could silently
        % disagree at the edge.
        threshPos = ea_sweetspot_exact_maxstat_threshold(maxRperm, alphaMaxstat, 'right');
        threshNeg = ea_sweetspot_exact_maxstat_threshold(minRperm, alphaMaxstat, 'left');

        thisMaxstatMap = nan(size(Remp));
        thisMaxstatMap(Remp >= threshPos) = Remp(Remp >= threshPos);
        thisMaxstatMap(Remp <= threshNeg) = Remp(Remp <= threshNeg);
        maxstatMap{side} = thisMaxstatMap;

        results.maxstat.threshPos{gi,side} = threshPos;
        results.maxstat.threshNeg{gi,side} = threshNeg;
        results.maxstat.thresholdedMap{gi,side} = thisMaxstatMap;

        [~, results.maxstat.pPos{gi,side}, results.maxstat.rankPos{gi,side}] = ea_sweetspot_nulldist_plot( ...
            maxRperm, maxRemp, 'Max R across all voxels', ...
            sprintf('Max-Statistic FWER (Positive, %s)', tag), 'right', useTailApprox);

        [~, results.maxstat.pNeg{gi,side}, results.maxstat.rankNeg{gi,side}] = ea_sweetspot_nulldist_plot( ...
            minRperm, minRemp, 'Min R across all voxels', ...
            sprintf('Max-Statistic FWER (Negative, %s)', tag), 'left', useTailApprox);

        fprintf('%s: voxelwise %d/%d pos + %d/%d neg voxels survive p<=%.3g | Eisenstein Q p=%.3g | max-stat pos p=%.3g, neg p=%.3g\n', ...
            tag, sum(voxSurvivePos), numel(voxSurvivePos), sum(voxSurviveNeg), numel(voxSurviveNeg), alphaVoxelwise, ...
            results.eisenstein.p{gi,side}, results.maxstat.pPos{gi,side}, results.maxstat.pNeg{gi,side});

        summaryLines{end+1} = sprintf('%s: voxelwise %d/%d pos + %d/%d neg voxels survive p<=%.3g | Eisenstein Q p=%.3g | max-stat pos p=%.3g, neg p=%.3g', ...
            tag, sum(voxSurvivePos), numel(voxSurvivePos), sum(voxSurviveNeg), numel(voxSurviveNeg), alphaVoxelwise, ...
            results.eisenstein.p{gi,side}, results.maxstat.pPos{gi,side}, results.maxstat.pNeg{gi,side}); %#ok<AGROW>
    end

    ea_sweetspot_vals2nii(train.space, voxThreshMap, outdir, sprintf('%s_voxelwise_p%s_group%d', basefname, alphaVoxelwiseStr, gi));
    ea_sweetspot_vals2nii(train.space, maxstatMap, outdir, sprintf('%s_maxstat_p%s_group%d', basefname, alphaMaxstatStr, gi));
end

readmeBody = sprintf([ ...
    'Post-hoc significance analysis of %s.\n\n', ...
    'Settings used for this run:\n', ...
    '  alphaVoxelwise  = %g\n', ...
    '  alphaEisenstein = %g\n', ...
    '  alphaMaxstat    = %g\n', ...
    '  useTailApprox   = %d\n\n', ...
    'Results:\n  %s\n\n', ...
    'Files added to this folder:\n', ...
    '  %s_voxelwise_p%s_group*.nii - per-voxel uncorrected permutation threshold, pos/neg\n', ...
    '  %s_maxstat_p%s_group*.nii   - max-statistic FWER-corrected threshold, pos/neg\n'], ...
    basefname, alphaVoxelwise, alphaEisenstein, alphaMaxstat, useTailApprox, ...
    strjoin(summaryLines, sprintf('\n  ')), basefname, alphaVoxelwiseStr, basefname, alphaMaxstatStr);
ea_sweetspot_readme_append(outdir, 'Post-hoc significance analysis (ea_sweetspot_permtest_stats)', readmeBody);
end

function thresh = ea_sweetspot_exact_maxstat_threshold(nullvals, alpha, tail)
% Exact-rank threshold matching the (count+1)/(Nperm+1) p-value formula used
% throughout this file and in ea_sweetspot_nulldist_plot -- the value at
% which a tied empirical statistic would first report p <= alpha, using
% only order statistics actually present in nullvals (unlike prctile, which
% interpolates between them).
%
% tail: 'right' (large nullvals are extreme, e.g. positive max-stat)
%       'left'  (small nullvals are extreme, e.g. negative max-stat)
%
% Returns NaN if no achievable rank gets p <= alpha with this many
% permutations (i.e. even the single most extreme null value isn't
% enough) -- comparisons against NaN are always false, so downstream
% thresholding correctly marks nothing as surviving.

Nperm = numel(nullvals);
switch tail
    case 'right'
        sorted = sort(nullvals, 'descend');
    case 'left'
        sorted = sort(nullvals, 'ascend');
end

pAtRank = ((1:Nperm) + 1) / (Nperm + 1); % p-value if the empirical stat ties the k-th most extreme null value
survivingIdx = find(pAtRank <= alpha, 1, 'last');
if isempty(survivingIdx)
    thresh = NaN;
else
    thresh = sorted(survivingIdx);
end
end
