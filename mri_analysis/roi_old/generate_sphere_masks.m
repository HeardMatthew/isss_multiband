%% generate_sphere_masks

[MNI_XYZ,Region] = xlsread([pwd filesep 'multi_roi_list_12subj.xls']);
Region = Region((~cellfun(@isempty, Region)));
sphere_radius = 2;

home = pwd;
cd /users/PAS1342/osu9912/fmri/isss_multiband_preprocessed_batch_03/data_14subjanalysis/ZG_03Nov17/design/multi/
Vmask = spm_vol('mask.nii');
[mask_mat, mask_xyz] = spm_read_vols(Vmask);
cd(home)

mask_idx = find(mask_mat);     %%% The indices of where mask=1
num_mask_voxels = length(mask_idx);
x_size = Vmask.dim(1);
y_size = Vmask.dim(2);
z_size = Vmask.dim(3);
x_coord_vec = 1:x_size;
y_coord_vec = 1:y_size;
z_coord_vec = 1:z_size;
%%%%%%%%%%%%%%% Now make a 3D mesh-grid of x,y and z coords
%%% x-coord is the row in this grid
x_coord_grid = x_coord_vec' * ones(1,y_size);
%%% y-coord is the col in this grid
y_coord_grid = ones(x_size,1) * y_coord_vec;
%%% Now stack these in the 3rd dimension, to make coords cubes
x_coord_cube = [];
y_coord_cube = [];
z_coord_cube = [];

for z_slice = 1:z_size
    x_coord_cube = cat(3,x_coord_cube,x_coord_grid);
    y_coord_cube = cat(3,y_coord_cube,y_coord_grid);
    z_coord_cube = cat(3,z_coord_cube,z_slice*ones(size(x_coord_grid)));
end

for ROI=1:numROI % For each ROI...
    this_vox = MNI_XYZ(ROI,:)';
    disp(Region{ROI})

    vox_idx = find(all(mask_xyz == this_vox));
    [center_x, center_y, center_z] = ind2sub(size(mask_mat), vox_idx);

    distance_cube = sqrt( (x_coord_cube - center_x).^2 + ...
        (y_coord_cube - center_y).^2 + ...
        (z_coord_cube - center_z).^2  );

    within_sphere = (distance_cube <= sphere_radius);
    within_sphere_and_mask = within_sphere.*mask_mat;
    
    Vroi = Vmask;
    Vroi.fname = [Region{ROI} '_sphere_.nii'];
    Vroi.private.dat.fname = [Region{ROI} 'sphere_.nii'];
    
    spm_write_vol(Vroi, within_sphere_and_mask);
    
end