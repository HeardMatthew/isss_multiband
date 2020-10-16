%% find_overlapping
% Creates list of overlapping and unique voxels in ROIs across two sets of
% ROIs. 

function find_overlapping(subj, study, rois, roi1, roi2)
%% check inputs
if length(subj) < 2
    error('Submit ALL subjects!')
end

if ~isstruct(subj) || ~isstruct(study)
    error('subj and study are both struct!')
end

if ~isnumeric(dd) || ~isnumeric(ss)
    error('dd and ss specify number of design and scan!')
end

if ~ismember(create, [0 1]) || ~ismember(plot, [0 1])
    error('flags to specify create and plot should be 0 or 1!')
end



%%% finding out how many voxels in vocode apsl appeared in gap PPV
clear all; close all; clc;


cd '/fMRI/vocode/OLD38/documents/summary/april2020/'

Vol_gap=spm_vol('gap_y35_frontal_R.nii');
Img_gap=spm_read_vols(Vol_gap);

Vol_vocode=spm_vol('vocode_old38_Frontal_R.nii');
Img_vocode=spm_read_vols(Vol_vocode);


num_x=size(Img_gap,1); num_y=size(Img_gap,2); num_z=size(Img_gap,3);

num_overlapping_voxels=0;
num_exclusiv_vocode=0;
num_exclusiv_gap=0;

gap_only=zeros(num_x, num_y, num_z); 
vocode_only=zeros(num_x, num_y, num_z); 
overlapping=zeros(num_x, num_y, num_z); 

for x=1:num_x
    for y=1:num_y
        for z=1:num_z
            
            if Img_vocode(x,y,z)==0 && Img_gap(x,y,z) ==1
                gap_only(x,y,z)=1; 
                num_exclusiv_gap=num_exclusiv_gap+1; 
                
            elseif  Img_vocode(x,y,z)==1 && Img_gap(x,y,z) ==0  
                vocode_only(x,y,z)=1; 
                num_exclusiv_vocode=num_exclusiv_vocode+1; 
                
            elseif  Img_vocode(x,y,z)==1 && Img_gap(x,y,z) ==1 
                overlapping(x,y,z)=1; 
                num_overlapping_voxels=num_overlapping_voxels+1; 
  
            end
        end
    end
    
    
end


V = Vol_gap; 

V.fname = 'gap_only_voxels.nii';  %%% Give this a new name
spm_write_vol(V,gap_only);


V.fname = 'vocode_only_voxels.nii';  %%% Give this a new name
spm_write_vol(V,vocode_only);


V.fname = 'vocode_gap_overlaping_voxels.nii';  %%% Give this a new name
spm_write_vol(V,overlapping);







end