function [MEP_Onset_time,MEP_Onset_index] = CON_Finder(EMG_wave,time,Threshold,direction,varargin)%,n,start,direction)
defaultn = 1;
p = inputParser;

validScalarPosNum = @(x) isnumeric(x);
addRequired(p,'EMG_wave',validScalarPosNum);
addRequired(p,'time',validScalarPosNum);
addRequired(p,'Threshold',validScalarPosNum); %MCD in this case
addRequired(p,'direction',@ischar);
addOptional(p,'n',defaultn,validScalarPosNum);
addParameter(p,'start',validScalarPosNum);
addParameter(p,'POI',validScalarPosNum);

parse(p,EMG_wave,time,Threshold,direction,varargin{:});
   
    
    if p.Results.direction == 'U'
        try
            for k = 1:(length(p.Results.EMG_wave))

                EMG_sec = [p.Results.EMG_wave(k:k+p.Results.n-1)]';

                Logical = EMG_sec > p.Results.Threshold; 
                Log_sum = sum(Logical);

                if Log_sum == p.Results.n
                   
                    MEP_Onset_index(k) = k;
                    MEP_Onset_time(k) = time(k);
                else
                    MEP_Onset_index(k) = nan;
                    MEP_Onset_time(k) = nan;

                end

            end
            
        catch
            for k = 1:(length(p.Results.EMG_wave)-p.Results.n)

                EMG_sec = [p.Results.EMG_wave(k:k+p.Results.n-1)]';

                Logical = EMG_sec > p.Results.Threshold; 
                Log_sum = sum(Logical);

                if Log_sum == p.Results.n
                   
                    MEP_Onset_index(k) = k;
                    MEP_Onset_time(k) = time(k);
                else
                    MEP_Onset_index(k) = nan;
                    MEP_Onset_time(k) = nan;

                end

            end
        end
            MEP_Onset_index = rmmissing(MEP_Onset_index');
            MEP_Onset_time = rmmissing(MEP_Onset_time');
%             ind = find(sum(MEP_Onset_index,2)==0) ;
%             MEP_Onset_index(ind,:)=[];
%             ind_2 = find(sum(MEP_Onset_time,2)==0);
%             MEP_Onset_time(ind_2,:)=[];

    else
        try
            for k = 1:(length(p.Results.EMG_wave))

                EMG_sec = [p.Results.EMG_wave(k:k+p.Results.n-1)]';

                Logical = EMG_sec < p.Results.Threshold; 
                Log_sum = sum(Logical);

                if Log_sum == p.Results.n
                   
                    MEP_Onset_index(k) = k;
                    MEP_Onset_time(k) = time(k);
                else
                    MEP_Onset_index(k) = nan;
                    MEP_Onset_time(k) = nan;

                end


            end
        catch
            for k = 1:(length(p.Results.EMG_wave)-p.Results.n)

                EMG_sec = [p.Results.EMG_wave(k:k+p.Results.n-1)]';

                Logical = EMG_sec < p.Results.Threshold; 
                Log_sum = sum(Logical);

                if Log_sum == p.Results.n
                   
                    MEP_Onset_index(k) = k;
                    MEP_Onset_time(k) = time(k);
                else
                    MEP_Onset_index(k) = nan;
                    MEP_Onset_time(k) = nan;

                end


            end

        end
        
    end

        MEP_Onset_index = rmmissing(MEP_Onset_index');
        MEP_Onset_time = rmmissing(MEP_Onset_time');

            if p.Results.POI
                MEP_Onset_index = p.Results.POI(MEP_Onset_index);
                MEP_Onset_time = p.Results.POI(MEP_Onset_time);
            end

if isempty(MEP_Onset_index)
    MEP_Onset_index = nan;
end

if isempty(MEP_Onset_time)
    MEP_Onset_time = nan;
end


end