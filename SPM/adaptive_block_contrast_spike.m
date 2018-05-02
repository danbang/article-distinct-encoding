function concatenate = adaptive_block_contrast_spike(v_contrasts,v_spike,n_block,n_motion,n_spike)
concatenate=[];
for j= 1:n_block
    if v_spike(j); 
        concatenate = [concatenate v_contrasts zeros(1,n_motion+n_spike)];
    else
        concatenate = [concatenate v_contrasts zeros(1,n_motion)]; 
    end
end
concatenate= [concatenate zeros(1,n_block)];
end