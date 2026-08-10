function ea_sweetspot_cap_parpool(maxWorkers)
% Ensures the current parallel pool (if any) has no more than maxWorkers
% workers before a parfor loop starts. Standalone-function counterpart to
% ea_sweetspot's capParpool method (which caps to obj.permtestMaxWorkers) --
% needed here because ea_sweetspot_permtest_tfce operates on a saved
% permtest file directly, without a live ea_sweetspot object to call a
% method on.

p = gcp('nocreate');
if isempty(p)
    parpool(maxWorkers);
elseif p.NumWorkers > maxWorkers
    delete(p);
    parpool(maxWorkers);
end
end
