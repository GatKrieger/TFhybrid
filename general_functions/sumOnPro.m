function sum_over_promoter = sumOnPro(checStruct, promoter_size,cerOpar)
global tss_struct
%%
if strcmp(cerOpar, 'cer')
    tss = tss_struct.cer;
elseif strcmp(cerOpar, 'par')
    tss = tss_struct.par; 
end


sample_names = fieldnames(checStruct.norm);

for s = 1 : length(sample_names)
    
    curr_sample = char(sample_names(s));
    curr_sample_data = checStruct.norm.(curr_sample);

    for i = 1:length(tss)

        curr_coor = tss(i,2:3);
        curr_chr = tss(i,1);

        if isnan(curr_chr)
            sop(i) = NaN;
            promoter_raw(i,1:promoter_size+1) =  NaN;

        elseif curr_coor(1) < curr_coor(2)   
            promoter_raw(i,1:promoter_size+1) = curr_sample_data{curr_chr}(curr_coor(1)- promoter_size: curr_coor(1));
            sop(i) = nansum(curr_sample_data{curr_chr}(curr_coor(1)- promoter_size: curr_coor(1)));
            
         elseif curr_coor(1) > curr_coor(2)  %reverse
            promoter_raw(i,1:promoter_size+1) = flip(curr_sample_data{curr_chr}(curr_coor(1):curr_coor(1)+promoter_size));
            sop(i) = nansum(curr_sample_data{curr_chr}(curr_coor(1):curr_coor(1)+promoter_size));


        end

    end
    sop(isnan(sop)) = 0;

    sum_over_promoter.(curr_sample) = sop;
    promoter_signal.(curr_sample) = promoter_raw;

end

%%
end
