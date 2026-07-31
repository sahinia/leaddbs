function results = ea_sweetspot_permtest_stats(permtestFile, alphaVoxelwise, alphaEisenstein, alphaMaxstat, useTailApprox)
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

train = load(permtestFile);
outdir = [fileparts(permtestFile), filesep];
[~, basefname] = fileparts(permtestFile);

results = struct();
results.sourcefile = permtestFile;
results.alphaVoxelwise = alphaVoxelwise;
results.alphaEisenstein = alphaEisenstein;
results.alphaMaxstat = alphaMaxstat;
results.useTailApprox = useTailApprox;

nsides = numel(train.space);

for gi = 1:size(train.Remp, 1)
    if all(cellfun(@isempty, train.Remp(gi,:)))
        continue
    end

    voxThreshMap = cell(1, nsides);
    maxstatMap = cell(1, nsides);

    for side = 1:nsides
        Remp = train.Remp{gi,side};   % Nvoxels x 1
        pemp = train.pemp{gi,side};   % Nvoxels x 1
        Rperm = train.Rperm{gi,side}; % Nperm x Nvoxels
        pperm = train.pperm{gi,side}; % Nperm x Nvoxels
        Nperm = train.Nperm;

        tag = sprintf('group%d side%d', gi, side);

        %% 1. Voxelwise uncorrected permutation threshold (per-voxel own null, pos/neg separate)
        RempRow = Remp'; % 1 x Nvoxels, for broadcasting against Rperm's rows
        countPos = sum(Rperm >= RempRow, 1);
        countNeg = sum(Rperm <= RempRow, 1);
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

        sigMaskPerm = pperm <= alphaEisenstein;
        logpPerm = -log10(max(pperm, eps));
        logpPerm(~sigMaskPerm) = 0;
        Qperm = ea_nansum(logpPerm, 2); % Nperm x 1

        results.eisenstein.Qemp{gi,side} = Qemp;
        results.eisenstein.Qperm{gi,side} = Qperm;

        [~, results.eisenstein.p{gi,side}, results.eisenstein.rank{gi,side}] = ea_sweetspot_nulldist_plot( ...
            Qperm, Qemp, 'Q = \Sigma -log_{10}(p)', ...
            sprintf('Eisenstein 2014 Omnibus Test (%s, \\alpha=%.3g)', tag, alphaEisenstein), 'right', useTailApprox);

        %% 3. Max-statistic FWER correction (pos/neg separate)
        % NOTE: ea_nanmax/ea_nanmin (ext_libs/nan) do NOT use MATLAB's own
        % max(A,[],dim) convention -- their 3-argument form is (a,dim,b) for
        % an ELEMENTWISE max/min of two same-sized arrays when dim is empty,
        % not "reduce along dim". The 2-argument form (a,dim) is what reduces
        % along a dimension.
        maxRperm = ea_nanmax(Rperm, 2); % Nperm x 1
        minRperm = ea_nanmin(Rperm, 2); % Nperm x 1
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
    end

    ea_sweetspot_vals2nii(train.space, voxThreshMap, outdir, sprintf('%s_voxelwise_p%.3g_group%d', basefname, alphaVoxelwise, gi));
    ea_sweetspot_vals2nii(train.space, maxstatMap, outdir, sprintf('%s_maxstat_p%.3g_group%d', basefname, alphaMaxstat, gi));
end
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
