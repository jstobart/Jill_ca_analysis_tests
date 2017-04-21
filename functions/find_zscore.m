function normTrace = find_zscore(trace, bgLevel, zIters, zSDs)

trace = trace-bgLevel;
raw_sig = trace;
idx = false(size(trace));

n = 0; 

while n < zIters
    
    sd_s = utils.nansuite.nanstd(trace(:));
    mean_s = utils.nansuite.nanmean(trace(:));
    idx = idx | (trace > (mean_s + zSDs*sd_s)) | ...
        (trace < (mean_s - zSDs*sd_s));
    
    raw_sig(idx) = NaN;
    
    n = n + 1; 
    
end

f0 = utils.nansuite.nanmean(raw_sig);
sigma = utils.nansuite.nanstd(raw_sig);

normTrace = (trace - f0) ./ sigma;

end