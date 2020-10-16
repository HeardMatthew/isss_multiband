%% combine_snr

function combine_snr(subj, study, dd, ss)

scan    = study.scan(ss); 
design  = study.design(dd); 
numsubj = length(subj); 

scanname = scan.runname(1:end-4); 
designname = [scanname '_' design.name]; 

for ns = 1:numsubj
    thissubj = subj(ns); 
    disp(thissubj.name)
    
    numruns  = thissubj.runs(ss); 
    
    dir_subj   = fullfile(study.path, 'data', thissubj.name); 
    dir_snr    = fullfile(dir_subj, 'SNR'); 
    
    dir_design = fullfile(dir_subj, 'design', designname); 
    mask = fullfile(dir_design, 'mask.nii'); 
    Vmask = spm_vol(mask); 
    ymask = spm_read_vols(Vmask); ymask = logical(ymask); 
    
    for nr = 1:numruns
        target = [study.prefix scan.runname num2str(nr)]; 
        filename = [target '_snr.nii']; 
        fullname = fullfile(dir_snr, target, 'SNR', filename); 
        
        V = spm_vol(fullname); 
        if ns == 1 && nr == 1
            data = nan([V.dim, numruns]);
        end
        
        thisdata = nan(V.dim); 
        
        temp = spm_read_vols(V); temp = temp(ymask); 
        thisdata(ymask) = temp; 
        data(:, :, :, nr) = thisdata; 
    end
    
    data_average = mean(data, 4);
    
    newname = [target(1:end-5) '_averaged_snr.nii'];
    V.fname = fullfile(dir_snr, newname); 
    spm_write_vol(V, data_average); 
end

end