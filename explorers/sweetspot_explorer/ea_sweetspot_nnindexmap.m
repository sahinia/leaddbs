function srcIdx = ea_sweetspot_nnindexmap(fromSpace, toSpace)
% Build a nearest-neighbor voxel correspondence between two nifti grids,
% without touching disk / SPM reslice. Mirrors what ea_conformspaceto with
% interp=0 (nearest neighbor, used throughout the sweetspot explorer)
% would produce geometrically, but computed once from header geometry
% alone so it can be reused to gather many maps' worth of values cheaply.
%
% fromSpace, toSpace - nii-like structs with .dim (1x3) and .mat (4x4),
%                       e.g. entries of obj.results.space.
%
% srcIdx - numel(toSpace grid) x 1. srcIdx(j) is the linear voxel index
%          into fromSpace's flattened image nearest to destination voxel
%          j, or NaN if voxel j falls outside fromSpace's grid entirely.

dimTo = toSpace.dim(1:3);
[X,Y,Z] = ndgrid(1:dimTo(1), 1:dimTo(2), 1:dimTo(3));
voxTo = [X(:)'; Y(:)'; Z(:)'; ones(1,numel(X))]; % 4 x Nvoxels, 1-based voxel coords

worldXYZ = toSpace.mat * voxTo; % voxel -> world mm

voxFrom = fromSpace.mat \ worldXYZ; % world mm -> source voxel (continuous)
voxFrom = round(voxFrom(1:3,:)); % nearest neighbor

dimFrom = fromSpace.dim(1:3);
valid = voxFrom(1,:)>=1 & voxFrom(1,:)<=dimFrom(1) & ...
        voxFrom(2,:)>=1 & voxFrom(2,:)<=dimFrom(2) & ...
        voxFrom(3,:)>=1 & voxFrom(3,:)<=dimFrom(3);

srcIdx = nan(numel(X),1);
srcIdx(valid) = sub2ind(dimFrom, voxFrom(1,valid), voxFrom(2,valid), voxFrom(3,valid));
