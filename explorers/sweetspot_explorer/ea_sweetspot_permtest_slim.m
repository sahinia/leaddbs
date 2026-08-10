function slimFile = ea_sweetspot_permtest_slim(permtestFile, slimFile)
% One-time conversion of a <ID>_permtest.mat file's Rperm/pperm (group 1
% only, matching predpermtest's existing restriction to train.Rperm{1,side})
% from Nperm x Nvoxels double arrays nested inside a {group,side} cell into
% flat, single-precision, top-level per-side variables in a new .mat file.
%
% Why this is necessary: matfile() only supports true partial (row-chunked)
% reads directly from disk for plain top-level array variables. Indexing
% into a cell array element (e.g. train.Rperm{1,side}) always materializes
% that element's FULL array in memory first -- there is no way to read a
% row subset of it lazily. At Nperm=5000 on a ~769K-voxel bilateral grid,
% Rperm+pperm together are ~60GB as double, which has crashed MATLAB
% outright on a 48GB-RAM machine via the blanket load() predpermtest used
% to do, even in a fresh session. Flattening Rperm/pperm out of the cell,
% one side at a time (so peak memory during conversion never holds more
% than one side's arrays at once), lets predpermtest read Nperm in batches
% afterwards instead of loading everything up front.
%
% Precision is also downcast double -> single here: correlation
% coefficients and p-values don't need double's ~16 significant digits;
% single's ~7 is far more than this analysis needs, and halves memory/disk
% size on top of the batching benefit.
%
% Skips the conversion (returns immediately) if slimFile already exists, so
% repeated predpermtest calls against the same training file (e.g. across
% different sigMode settings) only pay this cost once.
%
% Also skips it (returning permtestFile itself as slimFile) if permtestFile
% was already saved in this flat layout directly -- permtest() now does
% this by default for the common (single-group) case, making this whole
% conversion a no-op passthrough for any freshly-generated permtest file.
% This function, and the actual conversion below, still exist for older
% permtest.mat files saved before this change (or multi-group runs, which
% permtest() still saves in the original nested-cell format).

if ~exist('slimFile', 'var') || isempty(slimFile)
    [d, b] = fileparts(permtestFile);
    slimFile = fullfile(d, [b, '_slim.mat']);
end

vars = who('-file', permtestFile);
if ismember('Rperm_side1', vars)
    fprintf('%s is already in the flat, chunk-readable layout -- no conversion needed.\n', permtestFile);
    slimFile = permtestFile;
    return
end

if exist(slimFile, 'file')
    fprintf('Slim file already exists, reusing: %s\n', slimFile);
    return
end

fprintf('Creating slim (single-precision, chunk-readable) copy of %s...\n', permtestFile);

% Small variables (Nvoxels x 1, or smaller) -- safe to load whole
% regardless of Nperm/voxel count.
S = load(permtestFile, 'space', 'Remp', 'pemp', 'Nperm', 'corrtype', 'PermIdx', 'ID');
nsides = numel(S.space);

mIn = matfile(permtestFile);
mOut = matfile(slimFile, 'Writable', true);

mOut.space = S.space;
mOut.Remp = S.Remp;
mOut.pemp = S.pemp;
mOut.Nperm = S.Nperm;
mOut.corrtype = S.corrtype;
mOut.PermIdx = S.PermIdx;
mOut.ID = S.ID;
mOut.nsides = nsides;

for side = 1:nsides
    fprintf('  Converting side %d/%d (the memory-heavy step -- one side, one array, at a time)...\n', side, nsides);
    % mIn.Rperm(1,side) returns a 1x1 cell containing only this side's
    % array -- the sibling side never touches memory. R and p are also
    % handled one at a time (written and cleared before the other is
    % touched) rather than both held at once, since even a single side's
    % Rperm+pperm together can approach the physical RAM ceiling that
    % motivated this conversion in the first place.
    Rside = mIn.Rperm(1,side);
    mOut.(sprintf('Rperm_side%d', side)) = single(Rside{1});
    clear Rside
    pside = mIn.pperm(1,side);
    mOut.(sprintf('pperm_side%d', side)) = single(pside{1});
    clear pside
end

fprintf('Saved slim file to %s\n', slimFile);
end
