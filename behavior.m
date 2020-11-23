%%
% This code does the behavioral analyses described in the paper
% "Synchronization between keyboard typing and neural oscillations"
% It extracts RT and IKI data, prepares IKI data for Kernel Density
% Estimation analyses (KDE) and performs KDE.
% It outputs .mat files used for plotting and statistical analyses.
% This code applies to a matlab structure (named behavior in the following code) containing all the behavioral events
% occuring during the experiment. An example of such structure can be found
% at https://github.com/jduprez/Synchronization_btw_keyboard-oscillations


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

%% Initialize output matrices
% Create output matrix that will contain the average RT for each condition and outcome
nsub = 30;
ncond = 4; % word, pseudoword, sentence, pseudosentence
noutc = 3; % correct trials, corrected errors, other errors

RTmat = zeros(nsub, ncond, noutc);
RTcell = struct();

% Create output matrix for error details
errsubs = cell(nsub,1);
err_det = struct();

%% Loop through subjects

for subi = 1:nsub
    
    % Load behavioral data
    
    if subi<10
        
        load(strcat("P00",num2str(subi),"EEGTyping_behavior.mat"));
    else
        
        load(strcat("P0",num2str(subi),"EEGTyping_behavior.mat"));
    end
    
    %% Get RT data
    
    tempmat = zeros(max([behavior.Trial]), 4); % col1 = trial number, col2 = condition, col3 = outcome, col4 = RT
    % Conditions: 1 = word, 2 = pseudowords, 3 = sentences, 4 = pseudosentences
    % Outcome: 1 = correct trials, 2 = corrected errors, 3 = other errors
    
    for triali = min([behavior.Trial]):max([behavior.Trial])
        
        % Get rid of some buggy trials
        if subi == 26 && triali == 125
            continue
        end
        if subi == 29 && triali == 2
            continue
        end
        
        % switch conditions
        switch unique([behavior([behavior.Trial] == triali).Condition])
            
            case 1 % careful in the behavioral data, case 1 is free typing which is not analyzed in the paper
                
                continue
                
            case 2 % careful in the behavioral data, case 2 is pseudowords, becomes cond = 2
                
                if nanmean(unique([behavior([behavior.Trial] == triali).charmatch])) == 1 % correct trial
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 2;
                    tempmat(triali, 3) = 1;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                    
                elseif any([behavior([behavior.Trial] == triali).bs] == 1) % corrected error, bs is backspace
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 2;
                    tempmat(triali, 3) = 2;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                    
                else % other error
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 2;
                    tempmat(triali, 3) = 3;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                end
                
            case 3 % careful in the behavioral data, case 3 is pseudosentences, becomes cond = 4
                if nanmean(unique([behavior([behavior.Trial] == triali).charmatch])) == 1 % correct trial
                    disp('correct')
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 4;
                    tempmat(triali, 3) = 1;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                    
                elseif any([behavior([behavior.Trial] == triali).bs] == 1) % corrected error
                    disp('bs')
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 4;
                    tempmat(triali, 3) = 2;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                    
                else % other error
                    disp('err')
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 4;
                    tempmat(triali, 3) = 3;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                end
                
            case 4 % careful in the behavioral data, case 4 is words, becomes cond = 1
                if nanmean(unique([behavior([behavior.Trial] == triali).charmatch])) == 1 % correct trial
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 1;
                    tempmat(triali, 3) = 1;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                    
                elseif any([behavior([behavior.Trial] == triali).bs] == 1) % corrected error
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 1;
                    tempmat(triali, 3) = 2;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                    
                else % other error
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 1;
                    tempmat(triali, 3) = 3;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                end
                
            case 5 %% careful in the behavioral data, case 5 is sentences, becomes cond = 3
                if nanmean(unique([behavior([behavior.Trial] == triali).charmatch])) == 1 % correct trial
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 3;
                    tempmat(triali, 3) = 1;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                    
                elseif any([behavior([behavior.Trial] == triali).bs] == 1) % corrected error
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 3;
                    tempmat(triali, 3) = 2;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                    
                else % other error
                    tempmat(triali, 1) = triali;
                    tempmat(triali, 2) = 3;
                    tempmat(triali, 3) = 3;
                    templat = [behavior([behavior.Trial] == triali).latency];
                    tempmat(triali, 4) = templat(4) - templat(2);
                end
        end
    end
    
    RTcell(subi).allrt = tempmat; % store all RT of all subject in a structure
    
    for condi = 1:4
        for outi = 1:3
            temp = tempmat(tempmat(:,2) == condi & tempmat(:,3) == outi,:);
            avg = nanmean(temp, 1);
            RTmat(subi,condi, outi) = avg(4); % store subject averaged RT
        end
    end
    
    
    %% Compute condition-specific error-rates
    for triali = min([behavior.Trial]):max([behavior.Trial])
        
        err_det(triali).trial = triali; % store trial number
        err_det(triali).condition = unique([behavior([behavior.Trial] == triali).Condition]); % store condition
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 2 || unique([behavior([behavior.Trial] == triali).Condition]) == 4 % check if words/pseudowords
            
            if unique([behavior([behavior.Trial] == triali).charmatch]) == 1 % check if trial labelled correct (all characters matched)
                err_det(triali).acc = 1;
                err_det(triali).nBS = 0; % if trial correct, there's no backspace
            else
                err_det(triali).acc = 0;
            end
            err_det(triali).wordi = 1; % only one word since word/pseudoword condition
            err_det(triali).nerrWis = 0; % n error in a sentence if sentence
            err_det(triali).nBS = sum([behavior([behavior.Trial] == triali).bs]); % how many backspaces if any
            
            if sum([behavior([behavior.Trial] == triali).bs]) > 0
                err_det(triali).BSiki = mean([behavior([behavior.Trial] == triali & [behavior.bs] == 1).IKI]);
            else
                err_det(triali).BSiki = [];
            end
        else
            nerror = 0; % initialize n error at 0
            err_det(triali).nBS = sum([behavior([behavior.Trial] == triali).bs]); % n backspaces
            err_det(triali).BSiki = [behavior([behavior.Trial] == triali & [behavior.bs] == 1).IKI]; % backspaces IKI
            err_det(triali).wordi = max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount]);
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount]) % checks the number of erroneous words in sentences
                if unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 1
                    continue
                else
                    nerror = nerror+1;
                end
            end
            err_det(triali).nerrWis = nerror;
            if nerror>0
                err_det(triali).acc = 0;
            else
                err_det(triali).acc = 1;
            end
        end
    end
    
    % Put subject's results in a cell
    
    errsubs{subi} = err_det;
    
    
    %% Get IKI from behavior structure
    
    
    % Correct responses
    
    kpsw=1; kpss=1; knw=1; kns =1;
    
    
    for triali = min([behavior.Trial]):max([behavior.Trial]) % Use min instead of beginning at 1 because some trials were removed because of event problems
        
        % for pseudowords
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 2 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 1 ...
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4
            
            IKI_correct.pseudoW(kpsw, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            IKI_correct.pseudoW(kpsw, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            
            % create vector for fft
            latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
            vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
            vec4fft(latency_word - min(latency_word)+1) = 1;
            IKI_correct.pseudoW_vec4fft(kpsw, :) = {vec4fft};
            
            kpsw=kpsw+1;
            
            
        end
        
        % for normal words
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 4 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali ).charmatch]) == 1 ...
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4
            IKI_correct.W(knw, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            IKI_correct.W(knw, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            
            % create vector for fft
            latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
            vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
            vec4fft(latency_word - min(latency_word)+1) = 1;
            IKI_correct.W_vec4fft(knw, :) = {vec4fft};
            
            knw=knw+1;
        end
        
        % for pseudoword sentence
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 3
            
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if  unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 1 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                            
                    
                    IKI_correct.pseudoS(kpss, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                    IKI_correct.pseudoS(kpss, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                    
                    % create vector for fft
                    latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
                    vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
                    vec4fft(latency_word - min(latency_word)+1) = 1;
                    IKI_correct.pseudoS_vec4fft(kpss, :) = {vec4fft};
                    
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
                        
                        IKI_correct.S(kns, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        IKI_correct.S(kns, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        
                        % create vector for fft
                        latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
                        vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
                        vec4fft(latency_word - min(latency_word)+1) = 1;
                        IKI_correct.S_vec4fft(kns, :) = {vec4fft};
                        
                        kns = kns+1;
                    end
                end
            end
        end
        
        % disp(['correct_ ', num2str(triali)])
    end
    
    
    % Corrected errors (backspace)
    
    kpsw=1; kpss=1; knw=1; kns =1;
    
    
    for triali = min([behavior.Trial]):max([behavior.Trial])
        
        % check if backspace
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 2 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 0 ... % for pseudowords
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                && sum([behavior([behavior.Trial] == triali).bs]) > 0
            
            IKI_BSerr.pseudoW(kpsw, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            IKI_BSerr.pseudoW(kpsw, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            
            % create vector for fft
            latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
            vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
            vec4fft(latency_word - min(latency_word)+1) = 1;
            IKI_BSerr.pseudoW_vec4fft(kpsw, :) = {vec4fft};
            
            kpsw=kpsw+1;
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 4 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 0 ... % for pseudowords
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                && sum([behavior([behavior.Trial] == triali).bs]) > 0
            
            IKI_BSerr.W(knw, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            IKI_BSerr.W(knw, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            
            % create vector for fft
            latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
            vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
            vec4fft(latency_word - min(latency_word)+1) = 1;
            IKI_BSerr.W_vec4fft(knw, :) = {vec4fft};
            
            knw=knw+1;
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 3
            
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if  unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                            && sum([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).bs]) > 0
                        
                        IKI_BSerr.pseudoS(kpss, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        IKI_BSerr.pseudoS(kpss, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        
                        % create vector for fft
                        latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
                        vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
                        vec4fft(latency_word - min(latency_word)+1) = 1;
                        IKI_BSerr.pseudoS_vec4fft(kpss, :) = {vec4fft};
                        
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
                        IKI_BSerr.S(kns, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        IKI_BSerr.S(kns, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        
                        % create vector for fft
                        latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
                        vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
                        vec4fft(latency_word - min(latency_word)+1) = 1;
                        IKI_BSerr.S_vec4fft(kns, :) = {vec4fft};
                        
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
    
    
    for triali = min([behavior.Trial]):max([behavior.Trial])
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 2 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 0 ... % for pseudowords
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                && sum([behavior([behavior.Trial] == triali).bs]) == 0
            
            IKI_err.pseudoW(kpsw, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            IKI_err.pseudoW(kpsw, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            
            % create vector for fft
            latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
            vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
            vec4fft(latency_word - min(latency_word)+1) = 1;
            IKI_err.pseudoW_vec4fft(kpsw, :) = {vec4fft};
            
            kpsw=kpsw+1;
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 4 && length(unique([behavior([behavior.Trial] == triali).charmatch]) == 1)== 1 && unique([behavior([behavior.Trial] == triali).charmatch]) == 0 ... % for normal words
                && length([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                && sum([behavior([behavior.Trial] == triali).bs]) == 0
            
            IKI_err.W(knw, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            IKI_err.W(knw, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
            
            % create vector for fft
            latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
            vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
            vec4fft(latency_word - min(latency_word)+1) = 1;
            IKI_err.W_vec4fft(knw, :) = {vec4fft};
            
            knw=knw+1;
        end
        
        if unique([behavior([behavior.Trial] == triali).Condition]) == 3 % for pseudoword sentences
            
            for wordi = 1 : max([behavior([behavior.Trial] == triali & [behavior.wordcount] ~= 9999).wordcount])
                
                if isempty([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 % ensure charmatch is not empty
                    
                    if  unique([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).charmatch]) == 0 && length([behavior([behavior.Trial] == triali &[behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999 & ~isnan([behavior.CharDouble])  & [behavior.CharPos2] ~= 1).wordcount])>=4 ...
                            && sum([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi).bs]) == 0
                        
                        IKI_err.pseudoS(kpss, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        IKI_err.pseudoS(kpss, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        
                        % create vector for fft
                        latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
                        vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
                        vec4fft(latency_word - min(latency_word)+1) = 1;
                        IKI_err.pseudoS_vec4fft(kpss, :) = {vec4fft};
                        
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
                        IKI_err.S(kns, 1) = mean([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        IKI_err.S(kns, 2) = std([behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).IKI]);
                        % create vector for fft
                        latency_word = [behavior([behavior.Trial] == triali & [behavior.wordcount] == wordi & [behavior.wordcount] ~= 9999  & ~isnan([behavior.CharDouble]) & [behavior.CharPos2] ~= 1).latency];
                        vec4fft = zeros(1, max(latency_word)-min(latency_word)+1);
                        vec4fft(latency_word - min(latency_word)+1) = 1;
                        IKI_err.S_vec4fft(kns, :) = {vec4fft};
                        
                        kns = kns+1;
                    end
                end
            end
        end
        
        %disp(['other_errors ', num2str(triali)])
    end
    
    % Compute subject and condition averaged IKI for plotting and stats
    avgIKI_corr = zeros(30,4);
    avgIKI_BS = zeros(30,4);
    avgIKI_err = zeros(30,4);
    
    
    avgIKI_corr(subi, 1) = nanmean(IKI_correct.W(:,1), 1);
    avgIKI_corr(subi, 2) = nanmean(IKI_correct.pseudoW(:,1), 1);
    avgIKI_corr(subi, 3) = nanmean(IKI_correct.S(:,1), 1);
    avgIKI_corr(subi, 4) = nanmean(IKI_correct.pseudoS(:,1), 1);
    
    avgIKI_BS(subi, 1) = nanmean(IKI_BSerr.W(:,1), 1);
    avgIKI_BS(subi, 2) = nanmean(IKI_BSerr.pseudoW(:,1), 1);
    avgIKI_BS(subi, 3) = nanmean(IKI_BSerr.S(:,1), 1);
    avgIKI_BS(subi, 4) = nanmean(IKI_BSerr.pseudoS(:,1), 1);
    
    avgIKI_err(subi, 1) = nanmean(IKI_err.W(:,1), 1);
    avgIKI_err(subi, 2) = nanmean(IKI_err.pseudoW(:,1), 1);
    avgIKI_err(subi, 3) = nanmean(IKI_err.S(:,1), 1);
    avgIKI_err(subi, 4) = nanmean(IKI_err.pseudoS(:,1), 1);
    
    
    save([outfold 'P0' num2str(subi) '_IKI_FFT.mat'], 'IKI_correct', 'IKI_BSerr', 'IKI_err')
    
    %% Perform Kernel Density Estimation analyses
    
    fwhm = 1; % set full width at half maximum
    
    load(['P0' num2str(subi) '_IKI_FFT.mat'])
    results_kde_sub.correct(subi, 1) = {kde_letters_fun(IKI_correct.W, IKI_correct.W_vec4fft, fwhm)};
    results_kde_sub.correct(subi, 2) = {kde_letters_fun(IKI_correct.pseudoW, IKI_correct.pseudoW_vec4fft, fwhm)};
    results_kde_sub.correct(subi, 3) = {kde_letters_fun(IKI_correct.S, IKI_correct.S_vec4fft, fwhm)};
    results_kde_sub.correct(subi, 4) = {kde_letters_fun(IKI_correct.pseudoS, IKI_correct.pseudoS_vec4fft, fwhm)};
    results_kde_sub.BSerr(subi, 1) = {kde_letters_fun(IKI_BSerr.W, IKI_BSerr.W_vec4fft, fwhm)};
    results_kde_sub.BSerr(subi, 2) = {kde_letters_fun(IKI_BSerr.pseudoW, IKI_BSerr.pseudoW_vec4fft, fwhm)};
    results_kde_sub.BSerr(subi, 3) = {kde_letters_fun(IKI_BSerr.S, IKI_BSerr.S_vec4fft, fwhm)};
    results_kde_sub.BSerr(subi, 4) = {kde_letters_fun(IKI_BSerr.pseudoS, IKI_BSerr.pseudoS_vec4fft, fwhm)};
    results_kde_sub.err(subi, 1) = {kde_letters_fun(IKI_err.W, IKI_err.W_vec4fft, fwhm)};
    results_kde_sub.err(subi, 2) = {kde_letters_fun(IKI_err.pseudoW, IKI_err.pseudoW_vec4fft, fwhm)};
    results_kde_sub.err(subi, 3) = {kde_letters_fun(IKI_err.S, IKI_err.S_vec4fft, fwhm)};
    results_kde_sub.err(subi, 4) = {kde_letters_fun(IKI_err.pseudoS, IKI_err.pseudoS_vec4fft, fwhm)};
    
    % Prepare KDE data for plots and stats
    
    hz = linspace(0,100,1000); % Hz
    
    % Correct trials
    hz_correct = cellfun(@mean, results_kde_sub.correct,'UniformOutput',false); % average each cell
    % Extract words and compute mean.
    hz_correctW = cell2mat(hz_correct(:,1));
    hz_correctW = mean(hz_correctW, 1);
    % Extract pseudowords and compute mean
    hz_correctpseudoW = cell2mat(hz_correct(:,2));
    hz_correctpseudoW = mean(hz_correctpseudoW, 1);
    % Extract sentences and compute mean
    hz_correctS = cell2mat(hz_correct(:,3));
    hz_correctS = mean(hz_correctS, 1);
    % Extract pseudosentences and compute mean
    hz_correctpseudoS = cell2mat(hz_correct(:,4));
    hz_correctpseudoS = mean(hz_correctpseudoS, 1);
    
    % Corrected errors
    hz_BS = cellfun(@mean, results_kde_sub.BSerr,'UniformOutput',false); % average each cell
    hz_BSW = cell2mat(hz_BS(:,1));
    hz_BSW = mean(hz_BSW, 1);
    
    hz_BSpseudoW = cell2mat(hz_BS(:,2));
    hz_BSpseudoW = mean(hz_BSpseudoW, 1);
    
    hz_BSS = cell2mat(hz_BS(:,3));
    hz_BSS = mean(hz_BSS, 1);
    
    hz_BSpseudoS = cell2mat(hz_BS(:,4));
    hz_BSpseudoS = mean(hz_BSpseudoS, 1);
    
    % Errors
    hz_err = cellfun(@mean, results_kde_sub.err,'UniformOutput',false); % average each cell
    hz_errW = cell2mat(hz_err(:,1));
    hz_errW = mean(hz_errW, 1);
    
    hz_errpseudoW = cell2mat(hz_err(:,2));
    hz_errpseudoW = mean(hz_errpseudoW, 1);
    
    hz_errS = cell2mat(hz_err(:,3));
    hz_errS = mean(hz_errS, 1);
    
    hz_errpseudoS = cell2mat(hz_err(:,4));
    hz_errpseudoS = mean(hz_errpseudoS, 1);
    
    % Extract peak frequency
        
    % Correct trials
    hz_correct = results_kde_sub.correct;
    hz_correct = cellfun(@mean, hz_correct,'UniformOutput',false); % average each cell
    pfcorr = zeros(30,4);
    for subi = 1:nsub
        for condi = 1:4
            pfcorr(subi, condi) = hz(find(hz_correct{subi, condi}==max(hz_correct{subi, condi})));
            
        end
    end
    
    % BS errors
    hz_BS = results_kde_sub.BSerr;
    hz_BS = cellfun(@mean, hz_BS,'UniformOutput',false); % average each cell
    pfBS = zeros(30,4);
    for subi = 1:nsub
        for condi = 1:4
            pfBS(subi, condi) = hz(find(hz_BS{subi, condi}==max(hz_BS{subi, condi})));
            
        end
    end
    
    % Other errors
    
    hz_err = results_kde_sub.err;
    hz_err = cellfun(@mean, hz_err,'UniformOutput',false); % average each cell
    pferr = zeros(30,4);
    for subi = 1:nsub
        for condi = 1:4
            pferr(subi, condi) = hz(find(hz_err{subi, condi}==max(hz_err{subi, condi})));
            
        end
    end
    
    
    %%
    clear behavior err_det
    
    
    
    disp(strcat("Subject ", num2str(subi), " is done"))
end

% save
save([outfold, 'RTcell.mat'], 'RTcell')
save([outfold, 'RTmat.mat'], 'RTmat')
save([outfold, 'errsubs.mat'], 'errsubs')
save([outfold, 'avgIKI_corr.mat'], 'avgIKI_corr')
save([outfold, 'avgIKI_BS.mat'], 'avgIKI_BS')
save([outfold, 'avgIKI_err.mat'], 'avgIKI_err')
save([outfold, 'kderes4py.mat'], 'hz_correctW', 'hz_correctpseudoW', 'hz_correctS', 'hz_correctpseudoS', 'hz_BSW', 'hz_BSpseudoW', 'hz_BSS', 'hz_BSpseudoS', 'hz_errW' ,'hz_errpseudoW', 'hz_errS', 'hz_errpseudoS')
save([outfold, 'peakfreq_corr.mat'], 'pfcorr')
save([outfold, 'peakfreq_BS.mat'], 'pfBS')
save([outfold, 'peakfreq_err.mat'], 'pferr')


