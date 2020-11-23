%%
% This code does the EEG analyses described in the paper
% "Synchronization between keyboard typing and neural oscillations"
% It performs Generalized Eigen Decomposition (GED) and synchronization with
% typing (phase clustering).
% Note : GED requires inspection of the topography of the spatial filters
% to eliminate components showing artifact activity. Once this is done the
% best spatial filter is chosen based on the size of its associated
% eigenvalue to compute the component.
% This code applies to a EEGlab format files.

% Version 2 19/11/2020
% duprez.joan@gmail.com
% ## https://duprezjoan.wixsite.com/joanduprez ##

outfold = "/yourfold"; % insert output folder path

% Important :

% Condition 1 = Free typing task
% Condition 2 = PseudoWordList
% Condition 3 = PseudoWordSentenceList
% Condition 4 = NormalWordList
% Condition 5 = NormalWordSentenceList


nsub = 30;
ncond = 4; % word, pseudoword, sentence, pseudosentence
noutc = 3; % correct trials, corrected errors, other errors


% Initiate GED frequency parameters

min_frex_ged  =  1;
max_frex_ged  = 15;

% use linearly spaced frequencies
frex_ged = min_frex_ged:0.5:max_frex_ged;

% Initiate convolution frequency parameters

min_frex  =  1;
max_frex  = 50;
num_frex  = 60;
frex = logspace(log10(min_frex),log10(max_frex),num_frex); % logarithmically spaced frequency
nCycRange = [4 12];

% Initiate results matrices

% overall TF results stimulus onset epoching
tfged = zeros(size(frex_ged, 2),30,num_frex, 2201); % 29 ged frequencies, 30 participants*time points
% condition specific
tfged_cond = zeros(size(frex_ged, 2),30,num_frex, 2201, 5, 3); % 29 ged frequencies, 30 participants*time points* 4 conditions (1:W, 2:pW, 3:S, 4:pS) * 3 outcomes (Correct, corrected errors, other errors)


% overall ERP results
erpged = zeros(size(frex_ged, 2), 30, 2201); % 29 ged frequencies, 30 participants*time points
% condition specific erp results
erpged_cond = zeros(size(frex_ged, 2), 30, 2201, 5, 3); % 29 ged frequencies, 30 participants*time points* 4 conditions (1:W, 2:pW, 3:S, 4:pS) * 3 outcomes (Correct, corrected errors, other errors)

% matrix for significant clustering over all trials
sigfreq = zeros(30, 29);

%% First we perform GED and compute the spatial filters topographical map and save them in order to inspect offline and eliminate the bad filters
for subi = 1:nsub
    
    subnumber = subi;
    
    % Change subname according the number of the participants, important when calling for FT_Manual_Match
    
    if subi <= 9
        subname = ['P00' num2str(subi) 'EEGTyping_0.5_40_interp_rej'];
    else
        subname = ['P0' num2str(subi) 'EEGTyping_0.5_40_interp_rej'];
    end
    % Get ICs2remove
    
    IC2rm
    
    % Import EEG data (through EEGlab) --> EEGlab needs to be in the path
    % These EEG data are already preprocessed and the only step after
    % loading here is to remove ICA components previously defined and
    % stored in an ICs2remove .mat file
    %--------------------------------------------------------------------------
    
    eeglab
    
    % %     Import raw data
    EEG = pop_loadset([preprodatapath, subname '.set']);
    EEG = eeg_checkset( EEG );
    
    % remove ICA components
    EEG2 = pop_subcomp(EEG,ICs2remove{subnumber, 2});
    
    
    for freqi = 3:size(frex_ged,2)
        
        % perform GED on the continuous file
        % Needs the filterFGx function in the path to perform Narrow-band
        % filtering
        % Define S matrix
        
        tempdat   = EEG2.data; %
        
        td    = filterFGx2(tempdat,EEG.srate,frex(freqi),3); 
        td    = bsxfun(@minus,td,mean(td,2)); % subtract the mean (necessary)
        Smat = (td*td')/size(td,2);  % target frequency filtered channel covariance matrix
        
        % Define R matrix
        td = bsxfun(@minus, tempdat, mean(tempdat, 2));
        Rmat = (td*td')/size(td,2);  % Reference broadband channel covariance matrix
        
        % Regularization parameter
        
        shr = .005;
        Rmat = (1-shr)*Rmat + shr*mean(eig(Rmat))*eye(size(Rmat));
        
        % GED
        [original_evecsTh, original_evals] = eig(Smat,Rmat); 
        
        [evecsTh,evals] = eig(Smat,Rmat); % The columns of evecsTh are the spatial filters
        [~,sidx] = sort(diag(evals));
        evecsTh  = evecsTh(:,sidx); % Re-order evecsTh as a function of sidx which ordered the eigenvectors according to the size of their eigenvalue, does matter when the first vectors have the biggest eval !!
        mapsTh   = inv(evecsTh'); % forward model to see activation patterns
        
        % force sign of components --> so that topographical maps always show
        % positive values
        for ci=1:64
            [~,idx] = max(abs(mapsTh(:,ci))); % find strongest weight
            mapsTh(:,ci) = mapsTh(:,ci) * sign(mapsTh(idx,ci)); % force to positive sign
        end
        
        % Plot topographical map of the spatial filters, needs EEGlab as
        % well
        figure(freqi), clf
        
        for ploti = 1:15
            subplot(4,4, ploti)
            topoplot(mapsTh(:,end+1-ploti),EEG2.chanlocs,'numcontour', 6)
            title(num2str(64+1-ploti))
        end
        save ([outfold, num2str(subi) 'freq' num2str(freqi)], 'mapsTh', 'evecsTh', 'evals', 'original_evecsTh', 'original_evals') % keep all GED information !
        saveas(gcf,[outfold '/topomap_n' num2str(subi) 'freq' num2str(freqi) '.png'] ) % save topographical maps
        disp([ num2str(freqi) 'Hz done'])
    end
    disp([ 'subject ' num2str(subi) ' done'])
end

%% Before the next steps, the topography of the spatial filters needs to be inspected.
% When choice of the best spatial filter has to be forced (because the
% first component results from an artifact for example) the choice was
% stored in the GED_component_selection_05Hz_res.mat file as follows:
% if subi == 2
%     if freqgdi == 13
%         comp2keep = 63;
%     elseif freqgdi == 14
%         comp2keep = 63;
%     elseif freqgdi == 15
%         comp2keep = 63;
%     end
% end


%% Create subject*freq*comp matrix to plot average topomaps
avg_mapsTH = zeros(30, 29, 64);

for subi = 1:nsub
    for freqi = 1:29
        comp2keep = 64;
        load(['GEDresults_sub', num2str(subi), 'freq', num2str(freqi), '.mat'])
        GED_component_selection_05hz_res;
        avg_mapsTH(subi, freqi,: ) = squeeze(mapsTh(:, comp2keep));
    end
end


% Get the max eval for each sub frequency
% must put the 'topomapsGED_reg_0.5Hz' file in the path

mat4eval = zeros(30,29);

% frequency parameters for GED
min_frex  =  1;
max_frex  = 15;
num_frex  = 30;

frex = min_frex:0.5:max_frex;



for subi = [2:3, 5:18, 20:30]
    for freqi = 3:29
        load(strcat("GEDresults_sub",num2str(subi),"freq",num2str(freqi),".mat"))
        mat4eval(subi, freqi) = max(max(evals));
    end
end
save ([outfold, '/mat4eval.mat'], 'mat4eval')


%% Get the time indices needed for epoching and further eeg analyses

% Get the mean IKI per trial for each condition

% Loop over subjects
for subi = 1:nsub
    if subi <= 9
        filename = ['P00' num2str(subi) 'EEGTyping_behavior.mat'];
    else
        filename = ['P0' num2str(subi) 'EEGTyping_behavior.mat'];
    end
    load(filename)
    
    % Correct responses
    
    kpsw=1; kpss=1; knw=1; kns =1;
    
    if subi == 5
        firsttrial = 83;
    elseif  subi == 20
        firsttrial = 34;
    else
        firsttrial = 1;
    end
    
    for triali = firsttrial:max([behavior.Trial])
        
%         if subi == 3 && triali == 28 ||  subi == 3 && triali == 192% exception for that normal word trial that actually has 2 words
%             continue
%         end
%         
        % for pseudowords
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 2 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 1 ...
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4
            
            idx_correct.pseudoW(kpsw, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_correct.pseudoW(kpsw, 2) = {triali};
            
            kpsw=kpsw+1;
            
        end
        
        % for normal words
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 4 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali ).charmatch]) == 1 ...
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4
            
            idx_correct.W(knw, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_correct.W(knw, 2) = {triali};
            
            knw=knw+1;
        end
        
        % for pseudoword sentence
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 3
            
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 1 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                          
            idx_correct.pseudoS(kpss, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_correct.pseudoS(kpss, 2) = {triali};
                    
                    kpss = kpss+1;
                    end
                end
            end
        end
        
        % for normal sentences
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 5
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if  unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 1 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4
                        
            idx_correct.S(kns, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_correct.S(kns, 2) = {triali};
                        
                        kns = kns+1;
                    end
                end
            end
        end
        
        % disp(['correct_ ', num2str(triali)])
    end
    
    
    % Corrected errors (backspace)
    
    kpsw=1; kpss=1; knw=1; kns =1;
    
    
    for triali = firsttrial:max([behavior.Trial])
        
        % check if backspace
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 2 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 0 ... % for pseudowords
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                && sum([behavior([behavior.Trial] == triali).bs]) > 0
            
            idx_BSerr.pseudoW(kpsw, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_BSerr.pseudoW(kpsw, 2) = {triali};
            
            kpsw=kpsw+1;
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 4 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 0 ... % for pseudowords
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                && sum([behavior([behavior.Trial] == triali).bs]) > 0
            
            idx_BSerr.W(knw, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_BSerr.W(knw, 2) = {triali};
                        
            knw=knw+1;
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 3
            
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if  unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                            && sum([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).bs]) > 0
                        
            idx_BSerr.pseudoS(kpss, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_BSerr.pseudoS(kpss, 2) = {triali};
                        
                        kpss = kpss+1;
                    end
                end
            end
        end
        if unique([behavior([behavior.Trial] == triali).Condition]) == 5
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                            && sum([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).bs]) > 0
            
            idx_BSerr.S(kns, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_BSerr.S(kns, 2) = {triali};
                        
                        kns = kns+1;
                    end
                end
            end
        end
        
        %disp(['bacskpace_errors_', num2str(triali)])
    end
    
    
    % Other errors
    % The difference with the backspace errors is that in the if condition we
    % check that ---> sum([behavior([behavior.Trial] == triali).bs]) == 0
    
    kpsw=1; kpss=1; knw=1; kns =1;
    
    
    for triali = firsttrial:max([behavior.Trial])
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 2 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 0 ... % for pseudowords
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                && sum([behavior([behavior.Trial] == triali).bs]) == 0

            idx_err.pseudoW(kpsw, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_err.pseudoW(kpsw, 2) = {triali};
            
            kpsw=kpsw+1;
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 4 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 0 ... % for normal words
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                && sum([behavior([behavior.Trial] == triali).bs]) == 0
            
            idx_err.W(knw, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_err.W(knw, 2) = {triali};
                      
            knw=knw+1;
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 3 % for pseudoword sentences
            
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if  unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                            && sum([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).bs]) == 0
                        
            idx_err.pseudoS(kpss, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_err.pseudoS(kpss, 2) = {triali};
                        
                        kpss = kpss+1;
                    end
                end
            end
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 5 % for normal sentences
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                            && sum([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).bs]) == 0

            idx_err.S(kns, 1) = {[behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).urevent]};
            idx_err.S(kns, 2) = {triali};
                        
                        kns = kns+1;
                    end
                end
            end
        end
        
        %disp(['other_errors ', num2str(triali)])
    end
    
    save([outfold 'P0' num2str(subi) '_idx.mat'], 'idx_correct', 'idx_BSerr', 'idx_err')
    disp(['participant ' num2str(subi) ' done'])
end

%% Now we apply GED on the data and perform ERP analyses and time-frequency decomposition using complex Morlet wavelet convolution

for subi = 1:nsub
    
    if subi <= 9
        subname = ['P00' num2str(subi) 'EEGTyping_0.5_40_interp_rej'];
    else
        subname = ['P0' num2str(subi) 'EEGTyping_0.5_40_interp_rej'];
    end
    
    % Get ICs2remove
    
    IC2rm
    
    % Import EEG data y in the path (because in Preprocessing_script&fun)
    
    eeglab
    
    % Import raw data
    EEG = pop_loadset([preprodatapath, subname '.set']);
    EEG = eeg_checkset( EEG );
    
    % remove ICA components
    EEG2 = pop_subcomp(EEG,ICs2remove{subnumber, 2});
    
    
    for freqgdi = 1:size(frex_ged, 2)% (really long)
        
        comp2keep = 64; % default choice of spatial filter associated with the last component (the one with the highest eigenvalue)
        
        % force component selection in some cases
        GED_component_selection_05hz_res
        
        load(['GEDresults_sub' num2str(subnumber) 'freq' num2str(freqgdi) '.mat'])
        
        % create data
        GEDdata = (EEG2.data'*evecsTh(:,comp2keep))'; % gives a 1 component * n timepoints matrix
        
        
        %% load the preprocessed behavior data to get the  events
        load(['P0' num2str(subi), '_idx.mat'])
        if subi < 10
            load(['P00' num2str(subi), 'EEGTyping_behavior.mat'])
        else
            load(['P0' num2str(subi), 'EEGTyping_behavior.mat'])
        end
        
        %% Epoching around stimulus display for correct trials
        
        trial_idx = unique([idx_correct.W{:,2}, idx_correct.pseudoW{:,2},idx_correct.S{:,2},idx_correct.pseudoS{:,2}]); % get idx
        trial_idx=sort(trial_idx);
        trial_ev = zeros(size(trial_idx,2),2);
        for triali = 1:size(trial_idx,2)
            if isempty([behavior([behavior.Trial] == trial_idx(triali)).latency])
                continue
            end
            temp = [behavior([behavior.Trial] == trial_idx(triali)).latency];
            condition = [behavior([behavior.Trial] == trial_idx(triali)).Condition];
            trial_ev(triali,1) = temp(2);
            trial_ev(triali,2) = unique(condition);
        end
        trial_ev(any(trial_ev == 0, 2),:) = [];
        
        confcond = trial_ev(:,2); % create a condition vector
        
        tempdat = zeros(size(trial_ev,1), 2201);
        for triali =  1:size(trial_ev,1)
            if trial_ev(triali) < size(GEDdata,2) && trial_ev(triali)+1500< size(GEDdata,2)
                tempdat(triali,:) = GEDdata(trial_ev(triali)-700:trial_ev(triali)+1500);
            else
                break
            end
        end
        tempdat(any(tempdat == 0, 2),:) = [];
        if subnumber == 3
            tempdat(20,:) = [];
        end
        confcond = confcond(1:size(tempdat,1));
        
        tempdat=tempdat';
        
        %% Epoching around stimulus display for corrected errors trials
        
        trial_idx = unique([idx_BSerr.W{:,2}, idx_BSerr.pseudoW{:,2},idx_BSerr.S{:,2},idx_BSerr.pseudoS{:,2}]); % get idx
        trial_idx=sort(trial_idx);
        trial_ev = zeros(size(trial_idx,2),2);
        for triali = 1:size(trial_idx,2)
            if isempty([behavior([behavior.Trial] == trial_idx(triali)).latency])
                continue
            end
            temp = [behavior([behavior.Trial] == trial_idx(triali)).latency];
            condition_bserr = [behavior([behavior.Trial] == trial_idx(triali)).Condition];
            trial_ev(triali,1) = temp(2);
            trial_ev(triali,2) = unique(condition_bserr);
        end
        trial_ev(any(trial_ev == 0, 2),:) = [];
        
        confcond_bserr = trial_ev(:,2); % create a condition vector
        
        tempdat_bserr = zeros(size(trial_ev,1), 2201);
        for triali =  1:size(trial_ev,1)
            if trial_ev(triali) < size(GEDdata,2) && trial_ev(triali)+1500< size(GEDdata,2)
                tempdat_bserr(triali,:) = GEDdata(trial_ev(triali)-700:trial_ev(triali)+1500);
            else
                break
            end
        end
        tempdat_bserr(any(tempdat_bserr == 0, 2),:) = [];
        if subnumber == 3
            tempdat_bserr(20,:) = [];
        end
        confcond_bserr = confcond_bserr(1:size(tempdat_bserr,1));
        
        tempdat_bserr=tempdat_bserr';
        
        
        %% Epoching around stimulus display for other errors trials
        
        trial_idx = unique([idx_err.W{:,2}, idx_err.pseudoW{:,2},idx_err.S{:,2},idx_err.pseudoS{:,2}]); % get idx
        trial_idx=sort(trial_idx);
        trial_ev = zeros(size(trial_idx,2),2);
        for triali = 1:size(trial_idx,2)
            if isempty([behavior([behavior.Trial] == trial_idx(triali)).latency])
                continue
            end
            temp = [behavior([behavior.Trial] == trial_idx(triali)).latency];
            condition_err = [behavior([behavior.Trial] == trial_idx(triali)).Condition];
            trial_ev(triali,1) = temp(2);
            trial_ev(triali,2) = unique(condition_err);
        end
        trial_ev(any(trial_ev == 0, 2),:) = [];
        
        confcond_err = trial_ev(:,2); % create a condition vector
        
        tempdat_err = zeros(size(trial_ev,1), 2201);
        for triali =  1:size(trial_ev,1)
            if trial_ev(triali) < size(GEDdata,2) && trial_ev(triali)+1500< size(GEDdata,2)
                tempdat_err(triali,:) = GEDdata(trial_ev(triali)-700:trial_ev(triali)+1500);
            else
                break
            end
        end
        tempdat_err(any(tempdat_err == 0, 2),:) = [];
        if subnumber == 3
            tempdat_err(20,:) = [];
        end
        confcond_err = confcond_err(1:size(tempdat_err,1));
        
        tempdat_err=tempdat_err';
        
        %% Create a matrix containing all epochs and a corresponding condtion vector and outcome vector
        
        epochsdat = cat(2, tempdat, tempdat_bserr, tempdat_err);
        condition = cat (1, confcond, confcond_bserr, confcond_err);
        outcome = ([ones(1, size(tempdat, 2)), repmat(2, 1, size(tempdat_bserr, 2)), repmat(3, 1, size(tempdat_err, 2))])'; % 1 is correct trials, 2 is corrected errors, 3 is other errors
        
        %% ERP
        % Get overal ERP
        erpged(freqgdi, subi,:) = mean(epochsdat, 2);
        for outi = 1:3
            for condi = 1:5
                erpged_cond(freqgdi, subi, :, condi, outi) = mean(epochsdat(:, condition == condi & outcome == outi), 2);
            end
        end
        
        %% convolution
        
        basetime   = 200:500; % baseline period corresponds to -500 to - 200
        
        
        % wavelet parameters
        s = logspace(log10(nCycRange(1)),log10(nCycRange(2)),num_frex)./(2*pi.*frex); % standard deviation of the Gaussian (with a logarithmic increase in the number of cycles)
        t = -2:1/EEG.srate:2;
        halfwave = floor((length(t)-1)/2);
        
        % convolution parameters
        nData = size(epochsdat,2)*size(epochsdat,1);
        nWave = length(t);
        nConv = nData+nWave-1;
        
        % wavelets
        cmwX = zeros(num_frex,nConv);
        for fi=1:num_frex
            cmw = fft(  exp(1i*2*pi*frex(fi).*t) .* exp( (-t.^2)/(2*s(fi)^2) )  ,nConv);
            cmwX(fi,:) = cmw ./ max(cmw); % amplitude-normalize in the frequency domain
        end
        
        
        % Convolution and power extraction
        
        % fft of data, reshape in 1 time series of concatenated trials
        eegX = fft( reshape(epochsdat,1,nData) ,nConv); % 47 is the FCz electrode
        
        % loop over frequencies
        basebyfreq = zeros(num_frex, 1);
        for fi=1:num_frex
            
            % inverse fft of the multiplication of both signal and wavelets fourier power spectra
            as = ifft( eegX.*cmwX(fi,:) );
            
            % trim the result of convolution (since nConv = nData+nWave-1)
            as = as(halfwave+1:end-halfwave);
            
            % reshape in time*trials
            as = reshape(as,size(epochsdat, 1),size(epochsdat, 2));
            
            % condition-average baseline
            basepow = mean(mean( abs(as(basetime,:)).^2,2),1);
            basebyfreq(fi,:)= basepow; % store baseline for keystroke level TF
            % Mean power (dB) for conflict conditions
            tfged(freqgdi,subi,fi,:) = 10*log10( mean(abs(as).^2,2) ./ basepow );
            for outi = 1:3
                for condi = 1:5
                    tfged_cond(freqgdi,subi,fi,:, condi, outi) = 10*log10( mean(abs(as(:, condition == condi & outcome == outi)).^2,2) ./ basepow );
                end
            end
        end % end of frequency loop
        
        disp(['GED frequency ' num2str(frex_ged(freqgdi)) ' Hz of subject ' num2str(subi) ' is done'])
        
    end
    disp(['subject ' num2str(subi) ' is done'])
end
save([outfold 'tf_keylevel_results.mat'], 'erpged', 'erpgedcond', 'tfged', 'tfgedcond')


%% Now plot ERP and TF maps of 6 Hz and 8.5 Hz (change the indices if another frequency needed) and perform permutation analyses on the time-frequency data

% Need the topomaps to have been saved before ! Either using the topoplot
% function in that script, or after plotting with MNE in the python script.

% ERP
GEDfig = figure;

% 6 Hz
subplot(10,8,[1:3 9:11 17:19])
image(imread([outfold '/GEDtopomaps_withoutsensors_6.0Hz.png']) );
set(gca,'visible','off')

subplot(10,8, 25:27)
plot(-700:1500 ,squeeze(mean(erpged(11,:,:), 2)), 'Linewidth', 2)
ylabel('\muV', 'Fontsize', 16)
xlim([-400, 1200])
ylim([-0.06, 0.06])
set(gca, 'XTick', [-400, 0, 400, 800, 1200], 'YTick', [-0.05, 0, 0.05], 'Fontsize', 16)
xlabel('Time (ms)', 'Fontsize', 16)
line([822 822], [-0.06, 0.06],'Color','black','LineStyle','--', 'Linewidth', 1.5); % mean RT is 822 ms, so on the x axis it is 822 + 500

% 8.5 Hz
subplot(10,8,[5:7, 13:15, 21:23])
image(imread([outfold '/GEDtopomaps_withoutsensors_8.5Hz.png']) );
set(gca,'visible','off')

subplot(10,8, 29:31)
plot(-700:1500 ,squeeze(mean(erpged(16,:,:), 2)), 'Linewidth', 2)
ylabel('\muV', 'Fontsize', 16)
xlim([-400, 1200])
ylim([-0.06, 0.06])
set(gca, 'XTick', [-400, 0, 400, 800, 1200], 'YTick', [-0.05, 0, 0.05], 'Fontsize', 16)
xlabel('Time (ms)', 'Fontsize', 16)
line([822 822], [-0.06, 0.06],'Color','black','LineStyle','--', 'Linewidth', 1.5); % mean RT is 822 ms, so on the x axis it is 822 + 500


subplot(10,8,[41:43, 49:51, 57:59, 65:67, 73:75])

% Permutation analyses to show significant TF power changes from baseline

% For 6 Hz
times2save = 1:15:2200; % for temporal downsampling (at 67 Hz here) : otherwise the resulting matrices are way too big
tfged_downsampled = tfged(:,:,:,times2save);
nTimepoints = size(-700:15:1500, 2);
baseidx = dsearchn(times2save',[200 500]'); % corresponds to - 500 to - 200 ms prestimulus baseline

tfavg6 = squeeze(tfged_downsampled(11,:,:,:)); % take only the 6 Hz component

voxel_pval   = 0.05;
n_permutes = 500;

% initialize null hypothesis matrices
permuted_vals    = zeros(n_permutes,size(tfavg6, 1),num_frex,numel(times2save));



for permi=1:n_permutes
    
    for subi = 1:size(tfavg6, 1)
        cutpoint = randsample(2:nTimepoints-diff(baseidx)-2,1);
        permuted_vals(permi,subi,:,:) = tfavg6(subi,:,[cutpoint:end 1:cutpoint-1]);
    end
    disp(['permutation ', num2str(permi), ' is done']);
end
realmean = squeeze(mean(tfavg6, 1));
perm_mean = squeeze(mean(permuted_vals, 2));
zmap = (realmean-squeeze(mean(perm_mean,1)))./squeeze(std(perm_mean,1));
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;


% TF map 
climdb =[-0.25 0.25];
contourf(-700:15:1500,frex,realmean, 400,'linecolor','none')
set(gca, 'clim', climdb, 'FontSize',12,'yscale','log', 'xtick',[-400, 0, 400, 800, 1200],'ytick',logspace(log10(min(frex)),log10(max(frex)),6),'yticklabel',round(logspace(log10(min(frex)),log10(max(frex)),6)*10)/10, 'Fontsize', 16)
%title( '5 Hz GED - defined component power', 'FontSize', 20)
colormap jet
cbh = colorbar;
set(cbh, 'YTick', [-2.5, 2.5],'Fontsize', 13)
ylabel('Frequency (Hz)','Fontsize', 16)
xlabel('Time (ms)','Fontsize', 16)
xlim([-400, 1200])
cbh.Label.String = 'dB  from baseline';
cbh.FontSize = 18;
set(cbh, 'YTick', [-0.25, 0, 0.25])
hold on
line([822 822], [1 50],'Color','black','LineStyle','--', 'Linewidth', 1.5); % mean RT is 822 ms, so on the x axis it is 822 + 500
hold on
contour(-700:15:1500,frex,logical(zmapthresh), 1,'linecolor','w', 'linewidth', 2)


% Same for 8.5 Hz

%% Permutation analyses to show significant TF power changes from baseline
tfavg8_5 = squeeze(tfged_downsampled(16,:,:,:)); % take only the 8.5 Hz component

voxel_pval   = 0.05;
cluster_pval = 0.05;
n_permutes = 500;

% initialize null hypothesis matrices
permuted_vals    = zeros(n_permutes,size(tfavg8_5, 1),num_frex,numel(times2save));

for permi=1:n_permutes
    
    for subi = 1:size(tfavg8_5, 1)
        cutpoint = randsample(2:nTimepoints-diff(baseidx)-2,1);
        permuted_vals(permi,subi,:,:) = tfavg8_5(subi,:,[cutpoint:end 1:cutpoint-1]);
    end
    disp(['permutation ', num2str(permi), ' is done']);
end
realmean = squeeze(mean(tfavg8_5, 1));
perm_mean = squeeze(mean(permuted_vals, 2));
zmap = (realmean-squeeze(mean(perm_mean,1)))./squeeze(std(perm_mean,1));
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;

subplot(10,8,[45:47, 53:55, 61:63, 69:71, 77:79])

climdb =[-0.25 0.25];
contourf(-700:15:1500,frex,realmean, 400,'linecolor','none')
set(gca, 'clim', climdb, 'FontSize',12,'yscale','log', 'xtick',[-400, 0, 400, 800, 1200],'ytick',logspace(log10(min(frex)),log10(max(frex)),6),'yticklabel',round(logspace(log10(min(frex)),log10(max(frex)),6)*10)/10, 'Fontsize', 16)
%title( '5 Hz GED - defined component power', 'FontSize', 20)
colormap jet
cbh = colorbar;
set(cbh, 'YTick', [-2.5, 2.5],'Fontsize', 13)
ylabel('Frequency (Hz)','Fontsize', 16)
xlabel('Time (ms)','Fontsize', 16)
xlim([-400, 1200])
cbh.Label.String = 'dB  from baseline';
cbh.FontSize = 18;
set(cbh, 'YTick', [-0.25, 0, 0.25])
hold on
contour(-700:15:1500,frex,logical(zmapthresh), 1, 'linecolor','w', 'Linewidth', 2)
hold on
line([822 822], [1 50],'Color','black','LineStyle','--', 'Linewidth', 1.5); % mean RT is 822 ms, so on the x axis it is 822 + 500


filename = [outfold 'GED_ERP_TF_6&8.5Hz_perm.png'];

print(GEDfig, filename, '-dpng', '-r1000')
print(filename, '-dpng', '-r1000')


%% Compute clustering of phase at each frequency at the time of keystrokes
% Could be done inside the GED subject loop. However we perform permutation
% to zscore the clustering values at the same time, and this is very
% time-consuming.

nPerm = 500;
clustPerms = zeros(nPerm,1);
clustPermsZ = zeros(nPerm,1); % do the same but for ITPCz

for subi = 1:nsub
        
    if subi <= 9
        subname = ['P00' num2str(subi) 'EEGTyping_0.5_40_interp_rej'];
    else
        subname = ['P0' num2str(subi) 'EEGTyping_0.5_40_interp_rej'];
    end
        
    eeglab
    
    % Import raw data
    EEG = pop_loadset([preprodatapath, subname '.set']);
    EEG = eeg_checkset( EEG );
    % get ICs2remove
    IC2rm 
    
    % remove ICA components
    EEG2 = pop_subcomp(EEG,ICs2remove{subnumber, 2});
    
    % load indices for clustering
    load(['P0', num2str(subi), '_idx.mat'])
    
    for freqi = 1:29
        %% perform GED
                
        comp2keep = 64;
        
        % force component selection in some cases
        GED_component_selection_05hz_res
        
        load(['GEDresults_sub' num2str(subnumber) 'freq' num2str(freqi) '.mat'])
        
        % create data
        GEDdata = (EEG2.data'*evecsTh(:,comp2keep))'; % gives a 1 component * n timepoints matrix
        
        
        %% Hilbert transform
        % Filter-Hilbert
        hfiltdat = reshape(hilbert(filterFGx2(GEDdata,EEG.srate,frex(freqi),6)' ).', 1, EEG.pnts);
        
        %% take phase angle at each event
        phase_ang = angle(hfiltdat);
        
        %% Create the H0 distribution by permutation
        % MXC: permutation testing
        % get all time indices because we compute H0 on all the data and
        % not in a condition/outcome specific manner
        idx = [idx_correct.W{:,1},idx_correct.pseudoW{:,1}, idx_correct.S{:,1}, idx_correct.pseudoS{:,1}...
            idx_BSerr.W{:,1},idx_BSerr.pseudoW{:,1}, idx_BSerr.S{:,1}, idx_BSerr.pseudoS{:,1}...
            idx_err.W{:,1},idx_err.pseudoW{:,1}, idx_err.S{:,1}, idx_err.pseudoS{:,1}]; % get idx
        
        
        randcut = randperm(round(EEG.pnts*.9),nPerm) + round(EEG.pnts*.05);
        for permi=1:nPerm % cut-and-swap
            phase_angPerm = phase_ang([ randcut(permi):end 1:randcut(permi)-1 ]);
            clustPerms(permi) = abs(mean(exp(1i*(phase_angPerm(idx))))); % get clustering
            clustPermsZ(permi) = size(phase_angPerm,2)*(abs(mean(exp(1i*(phase_angPerm(idx)))))^2);
        end
        
        %% for correct trials
        
        % for words
        idx = [idx_correct.W{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.correct(freqi).W(subi,1) = clust;
        clustering.correct(freqi).W(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.correctZ(freqi).W(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.correctZ_z(freqi).W(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        
        % for pseudowords
        idx = [idx_correct.pseudoW{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.correct(freqi).pseudoW(subi,1) = clust;
        clustering.correct(freqi).pseudoW(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.correctZ(freqi).pseudoW(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.correctZ_z(freqi).pseudoW(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        % for sentences
        idx = [idx_correct.S{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.correct(freqi).S(subi,1) = clust;
        clustering.correct(freqi).S(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.correctZ(freqi).S(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.correctZ_z(freqi).S(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        % for pseudo sentences
        
        idx = [idx_correct.pseudoS{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.correct(freqi).pseudoS(subi,1) = clust;
        clustering.correct(freqi).pseudoS(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.correctZ(freqi).pseudoS(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.correctZ_z(freqi).pseudoS(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        %% for BS errors
        
        % for words
        idx = [idx_BSerr.W{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.BSerr(freqi).W(subi,1) = clust;
        clustering.BSerr(freqi).W(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.BSerrZ(freqi).W(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.BSerrZ_z(freqi).W(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        % for pseudowords
        idx = [idx_BSerr.pseudoW{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.BSerr(freqi).pseudoW(subi,1) = clust;
        clustering.BSerr(freqi).pseudoW(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.BSerrZ(freqi).pseudoW(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.BSerrZ_z(freqi).pseudoW(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        % for sentences
        idx = [idx_BSerr.S{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.BSerr(freqi).S(subi,1) = clust;
        clustering.BSerr(freqi).S(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.BSerrZ(freqi).S(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.BSerrZ_z(freqi).S(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        % for pseudo sentences
        
        idx = [idx_BSerr.pseudoS{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.BSerr(freqi).pseudoS(subi,1) = clust;
        clustering.BSerr(freqi).pseudoS(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.BSerrZ(freqi).pseudoS(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.BSerrZ_z(freqi).pseudoS(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        %% for other errors
        
        % for words
        idx = [idx_err.W{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.err(freqi).W(subi,1) = clust;
        clustering.err(freqi).W(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.errZ(freqi).W(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.errZ_z(freqi).W(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        % for pseudowords
        idx = [idx_err.pseudoW{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.err(freqi).pseudoW(subi,1) = clust;
        clustering.err(freqi).pseudoW(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.errZ(freqi).pseudoW(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.errZ_z(freqi).pseudoW(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        % for sentences
        idx = [idx_err.S{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.err(freqi).S(subi,1) = clust;
        clustering.err(freqi).S(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.errZ(freqi).S(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.errZ_z(freqi).S(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        % for pseudo sentences
        
        idx = [idx_err.pseudoS{:,1}]; % get idx
        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        zclust = size(idx, 2) * (clust^2); % corrected clustering with Rayleigh's Z
        clustering.err(freqi).pseudoS(subi,1) = clust;
        clustering.err(freqi).pseudoS(subi,2) = zclust;
        
        % convert to z-score distance to H0 distribution
        clustering.errZ(freqi).pseudoS(subi,1) = (clust-mean(clustPerms)) / std(clustPerms);
        clustering.errZ_z(freqi).pseudoS(subi,1) = (clust-mean(clustPermsZ)) / std(clustPermsZ);
        
        
        % Get overall clustering significance
        
              
        % for words
        idx = [idx_correct.W{:,1}, idx_correct.pseudoW{:,1}, idx_correct.S{:,1}, idx_correct.pseudoS{:,1}...
            idx_BSerr.W{:,1}, idx_BSerr.pseudoW{:,1}, idx_BSerr.S{:,1}, idx_BSerr.pseudoS{:,1}...
            idx_err.W{:,1}, idx_err.pseudoW{:,1}, idx_err.S{:,1}, idx_err.pseudoS{:,1}]; % get idx
        
        % Create the H0 distribution by permutation
        randcut = randperm(round(EEG.pnts*.9),nPerm) + round(EEG.pnts*.05);
        for permi=1:nPerm % cut-and-swap
            phase_angPerm = phase_ang([ randcut(permi):end 1:randcut(permi)-1 ]);
            clustPerms(permi) = abs(mean(exp(1i*(phase_angPerm(idx))))); % get clustering
        end

        phase_ang_idx = phase_ang(idx); % get phase angles
        phase_clust = exp(1i*(phase_ang_idx)); % put in the complex plane
        clust = abs(mean(phase_clust)); % get clustering
        sigfreq(subi, freqgdi) = exp(-size(idx,2)*(clust^2));% get significance
        
        
        disp([num2str(frex(freqi)), ' Hz of subject ', num2str(subi), ' is done'])
    end
end

% prepare clustering data for plotting and stats

clust_mat = zeros(3,29,4); % 3 outcome (corr, bs, err) * 29 freq * 4 cond (W, pW, S, pS), * 2 measures (ITPC, ITPCperm)

for subi = 1:nsub
    for freqi = 3:29
        
        clust_mat(1, freqi, 1) = mean(nonzeros(clustering.correct(freqi).W(:,1))) ;
        clust_mat(1, freqi, 2) = mean(nonzeros(clustering.correct(freqi).pseudoW(:,1))) ;
        clust_mat(1, freqi, 3) = mean(nonzeros(clustering.correct(freqi).S(:,1))) ;
        clust_mat(1, freqi, 4) = mean(nonzeros(clustering.correct(freqi).pseudoS(:,1))) ;
        
        clust_mat(2, freqi, 1) = mean(nonzeros(clustering.BSerr(freqi).W(:,1))) ;
        clust_mat(2, freqi, 2) = mean(nonzeros(clustering.BSerr(freqi).pseudoW(:,1))) ;
        clust_mat(2, freqi, 3) = mean(nonzeros(clustering.BSerr(freqi).S(:,1))) ;
        clust_mat(2, freqi, 4) = mean(nonzeros(clustering.BSerr(freqi).pseudoS(:,1))) ;
        
        clust_mat(3, freqi, 1) = mean(nonzeros(clustering.err(freqi).W(:,1))) ;
        clust_mat(3, freqi, 2) = mean(nonzeros(clustering.err(freqi).pseudoW(:,1))) ;
        clust_mat(3, freqi, 3) = mean(nonzeros(clustering.err(freqi).S(:,1))) ;
        clust_mat(3, freqi, 4) = mean(nonzeros(clustering.err(freqi).pseudoS(:,1))) ;
        
        
    end
end

save([outfold 'clust_mat.mat'], 'clust_mat')
save([outfold 'sigfreq.mat'], 'sigfreq')


