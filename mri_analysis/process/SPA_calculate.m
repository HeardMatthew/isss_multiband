%% SPA_calculate
% Creates Sum-Positive-Area underneath the curve for second-level FIR
% analysis. 
% MM/DD/YY -- CHANGELOG
% 04/07/20 -- Initialized changelog. Cloned from original hybrid_isss batch
%   for new analysis. Updated format to match new universal strategy. 

function SPA_calculate(varargin)
%% Check input
subj = varargin{1}; 
study = varargin{2}; 
dd = varargin{3}; 
ss     = varargin{4}; 

if ~isstruct(subj) || length(subj) ~= 1
    error('Input ("subj") where subj is a SINGLE struct')
end

if ~isstruct(study)
    error('Input ("study") which has experiment info')
end

if ~isnumeric(dd)
    error('Input ("dd") which specifies which designs (#) should be used')
end

if ~isnumeric(ss)
    error('Input ("ss") which specifies scan!')
end

%% Pathing
cd ..
dir_batch = pwd; 
dir_subj  = fullfile(study.path, 'data', subj.name); 
dir_design_root = fullfile(dir_subj, 'design');

%% Parameters
scan = study.scan(ss); 
scanname = scan.runname(1:end-4); 
thisdesign = study.design(dd); 
dir_thisdesign = fullfile(dir_design_root, [scanname '_' thisdesign.name]); 

numruns = subj.runs(ss);     
numcons = length(thisdesign.cond); 
numepis = 8/scan.TR; 

%% Load necessary data
disp(['Making contrast images for subject ' subj.name]);
disp(['Working on ' scanname '...'])
    
spmmat = fullfile(dir_thisdesign, 'SPM.mat');
load(spmmat)
    
mask = fullfile(dir_thisdesign, 'mask.nii');
Vmask = spm_vol(mask); % loads header for mask
mask_matrix = spm_read_vols(Vmask);
mask_inds = find(mask_matrix);     %%% The indices of where mask=1
num_mask_voxels = length(mask_inds);
    
[x_in_mask,y_in_mask,z_in_mask] = ind2sub(size(mask_matrix),mask_inds);
XYZ = [x_in_mask,y_in_mask,z_in_mask]';
    
for rr = 1:numruns % For each run...
    runtag = ['Sn(' num2str(rr) ')']; 

    for cc = 1:numcons
        cond_name = thisdesign.cond{cc}; 
            
        thiscon = contains(SPM.xX.name, [runtag ' ' cond_name]);
        thesebetas = cell(length(find(thiscon)), 1);
            
        if ~isempty(thesebetas)
            idx = 1;
            for vb = find(thiscon)
                thesebetas{idx} = fullfile(dir_thisdesign, SPM.Vbeta(vb).fname);
                idx = idx + 1;
            end

            betas_all = spm_get_data(thesebetas,XYZ);

            %% calculating area
            AUE_for_this_cond = zeros(1,num_mask_voxels);
            for voxel = 1:num_mask_voxels
                area_all = 0;
                for TR = 1:numepis-1 % minus one because last sample calls on TR end-1 and end
                    y1=betas_all(TR,voxel); y2=betas_all(TR+1,voxel);

                    % If both are above zero...
                    if y1 > 0 && y2 > 0 % trapezoid shape
                        area= (y1+ y2)/2;

                    % If both are below zero...
                    elseif y1 < 0 && y2 < 0
                        % area= -1*(y1+ y2)/2;
                        area=0;

                    % If 1 is above and 2 is below
                    elseif y1 > 0 && y2 < 0
                        b=y1;
                        a=y2-y1;
                        x=-b/a;
                        area=(x*y1/2);%+ ((1-x)*y2/2);

                    elseif y1 < 0 && y2 > 0
                        b=y1;
                        a=y2-y1;
                        x=-b/a;
                        %  area=(x*y1/2)+ ((1-x)*y2/2);
                        area=((1-x)*y2/2);

                    end

                    area_all=area_all+area;
                end

                AUE_for_this_cond(1,voxel)=area_all;

            end

            Vol = fullfile(dir_thisdesign, 'beta_0001.nii'); Vol = spm_vol(Vol);
            Img = spm_read_vols(Vol);
            Vol.fname = fullfile(dir_thisdesign, ['AUE_' cond_name '_run' num2str(rr) '.nii']);

            for voxel=1:num_mask_voxels
                Img(mask_inds(voxel))=AUE_for_this_cond(1,voxel);
            end

            spm_write_vol(Vol, Img);

        end % if empty line

    end
    
end

end
