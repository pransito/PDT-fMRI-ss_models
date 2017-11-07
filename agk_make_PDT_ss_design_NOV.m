% agk_make_PDT_ss_design_ed_03
% for one subject
% makes the ss design of the PDT task based on the Tom et al. (2007) paper
% pictures here are put in as categorical vectors
% to allow a cue reactivity contrast (between categories modeling)

% INPUT
% cur_sub       : name of current subject
% which_folders : name of folder with niftis of task (pciks first) ???
% aggr          : aggregation level gambling matrix (we usually take 3)
% run_it        : should the job run? (will be only saved if not)
% acc_rec       : should accept be included as a factor in model
% expl_mask     : usually use a gray matter mask here
% cur_tmpl      : the ss job template to be used here, cell, first field is
%                 the one without accept reject as factor and second with
% ow            : overwrite exisiting results files? (0,1)
% physio_inc    : should physio correction be incorp. (diff res folder)

% OUTPUT
% error_message: success message or explanation what hasn't worked

% IMPLICIT OUTPUT
% writes the ss model into the results_ssdir in the subject folder

% THIS MODEL ('PDT_ss_design_edfix_03')
% 4 task on: PIC, PIC+GAM, PIC+GAM+OPT
% 4 param modulators in PIC: category dummy coded (neu, gam, neg, pos)
% 3 param mod in PIC+GAM: gain, loss, ed & optional acceptance
% no param mod for PIC+GAM+OPT

% DETAILS
% preproc_job is needed to get the slice timing parameters
% (to create time bins), will be looked for in the nifti folder provided
% will look for "preprocessing..." pattern
% abs values of losses are used
% nifti data need to be SPM12 preprocessed
% takes in an aggregation factor: e.g. 1 means no change; 2 means, the 12
% by 12 matrix will be aggregated to a 6-by-6 matrix, and so on
% the picture category is modeled by three dummy variables (i.e. all 0 is
% neutral)
% EDABS still to aggregate (6 levels) and mean center

% AUTHORSHIP
% author: Alexander Genauck
% revision: Katharina Dücker
% date  : 23.01.2017, 27.10.2017, 2.11.2017
% email : alexander.genauck@charite.de

function [error_message] = agk_make_PDT_ss_design_NOV(cur_sub,cur_tmpl,aggr, ...
    run_it,acc_rec,expl_mask,ow,physio_inc,tr)
try
    %% PREPARATIONS
    % change into subs directory
    root_dir = pwd;
    cd(cur_sub)
    
    % name of this analysis determines the results directory
    if physio_inc
        results_ssdir = 'PDT_ss_design_NOV_phys';
    else
        results_ssdir = 'PDT_ss_design_NOV';
    end
    
    % messaging
    disp(['Estimating subject ' cur_sub ' using ' results_ssdir])
    
    % gain and loss ranges
    gain_min = 14;
    gain_max = 36;
    loss_max = 18;
    loss_min = 7;
    
    % for euclidean distance
    vec = [2;1;0];   % slope vector [gain; loss; 0]
    sp  = [26;13;0]; % sp support point; point on the diagonal [gain; loss; 0]
    
    %% FILLING BATCH
    % preparing the design information variables
    names       = {'Pic.on','Pic.gam.on','Pic.gam.opt.on'};
    
    % structures for onset and duration with subfields pic only, pic
    % and gamble and options
    % prepare cells for onsets and durations
    onsets          = cell(1,length(names));
    durations       = cell(1,length(names));
    
    pmod      = struct('name',{},'poly',{},'param',{});
    pic_pres  = {'gam', 'neg', 'pos'}; % pictures presented
    pic_gam   = {'gain', 'loss', 'ed'}; % picture + gamble on
    
    % if accepted included in model
    if acc_rec
        pic_gam = [pic_gam, {'acc'}];
        pic_gam = [pic_gam, {'accXgam'}];
        pic_gam = [pic_gam, {'accXneg'}];
        pic_gam = [pic_gam, {'accXpos'}];
    end
    
    % parametric modulators
    % 1. category while picture presented
    for pp = 1:length(pic_pres)
        pmod(1).name{pp} = pic_pres{pp};
        pmod(1).poly{pp} = 1;
    end
    pmod(1).param = cell(1,length(pic_pres));
    
    % 2. gain, loss, ed of current gamble
    for pp = 1:length(pic_gam)
        pmod(2).name{pp} = pic_gam{pp};
        pmod(2).poly{pp} = 1;
    end
    pmod(2).param = cell(1,length(pic_gam));
    
    % 3. no params for option on
    pmod(3).name = {};
    pmod(3).poly = {};
    pmod(3).poly = {};
    
    % get behav data
    behav_data = 0;
    try
        cd('Behav\PDT')
        
        if strcmp(cur_sub,'VPPG0115')
            % two p struct folders
            p_structs = cellstr(ls('0115*'));
            % go into the second folder
            cd(p_structs{2})
        end
        
        load(ls('P*'))
        behav_data = 1;
    catch
        behav_data = 0;
        error_message = [results_ssdir ' ' cur_sub ' No P struct found.'];
        cd(root_dir)
        return
    end
    
    % calculate the aggregated possible values
    cur_gains = (gain_min:2:gain_max);
    cur_losss = loss_min:loss_max    ;
    [vg,osg,nsg]  = agk_downsample_steps(cur_gains,aggr);
    gain_min  = min(vg);
    gain_max  = max(vg);
    [vl,osl,nsl]  = agk_downsample_steps(cur_losss,aggr);
    loss_min  = min(vl);
    loss_max  = max(vl);
    
    %getting the gain mean and the loss mean
    if strcmp(cur_sub,'VPPG0115')
        % getting the gain mean and the loss mean
        % special case: VPPG0115
        
        % gain and loss of second run
        start_run2 = length(P.cur.choice) - 100;
        % cut cur.choice and gamble (delete part that belongs to first
        % run)
        P.cur.choice           = P.cur.choice(start_run2:end);
        P.cur.gamble           = P.cur.gamble(start_run2:end);
        
        all_gains2 = [];
        all_losss2 = [];
        for r2 = 1:length(P.cur.choice)
            all_gains2(r2) = str2double(cell2mat(P.gain.strings(P.cur.gamble{r2}(1))));
            all_losss2(r2) = str2double(cell2mat(P.loss.strings(P.cur.gamble{r2}(2))));
        end
        % load p_struct of the first run
        % gain and loss of first run
        cd(root_dir)
        cd(cur_sub)
        cd('Behav\PDT')
        cd(p_structs{1})
        load(ls('P*'))
        P.cur.choice = P.cur.choice(1:101);
        P.cur.gamble = P.cur.gamble(1:101);
        
        P.cur_iti_dur_corr = P.cur_iti_dur_corr(1:101);
        
        for r1 = 1:length(P.cur.choice)
            all_gains1(r1) = str2double(cell2mat(P.gain.strings(P.cur.gamble{r1}(1))));
            all_losss1(r1) = str2double(cell2mat(P.loss.strings(P.cur.gamble{r1}(2))));
        end
        
        % concatenate both runs
        all_gains = [all_gains1 all_gains2];
        all_losss = [all_losss1 all_losss2];
        
    else
        
        for hh = 1:length(P.cur.choice)
            all_gains(hh) = str2double(cell2mat(P.gain.strings(P.cur.gamble{hh}(1))));
            all_losss(hh) = str2double(cell2mat(P.loss.strings(P.cur.gamble{hh}(2))));
        end
        
    end
    mean_gain = mean(all_gains);
    mean_loss = mean(all_losss);
    
    
    tworuns   = isfield(P.t,'triggerscannertpostpause');
    
    % check if this is two runs and get onset correction
    % added: length of each run + fill in pmods
    if tworuns
        % new starting point in time for second run
        time_corr_second_run_fixed = (P.t.triggerscannert0 - ...
            P.t.triggerscannertpostpause);
        % length of first run is always 101
        length_run = 101;
        num_run    = 2;
    else
        time_corr_second_run_fixed = 0;
        length_run                 = length(P.cur.choice);
        num_run                    = 1;
    end
    
    % prepare cell templates to fill in onsets, durations and param mod
    onsets_cell    = cell(1,num_run);     % onset cell -> pic, p+g, opt
    durations_cell = cell(1,num_run);     % durations cell
    pmod_cell      = cell(1,num_run);
    
    for nn = 1:num_run
        onsets_cell{nn}    = onsets;
        durations_cell{nn} = durations;
        pmod_cell{nn}      = pmod;
        % cell array for picture on only: matrix with zeros (???)
    end
    
    % some cleaning (only needed when building this)
    clear missing
    clear trial
    
    %% FILL ONSETS, PMODS
    % try
    length_choices = length(P.cur.choice);
    
    if strcmp(cur_sub,'VPPG0115')
        cd(root_dir)
        cd(cur_sub)
        cd('Behav\PDT')
        cd(p_structs{2})
        e = load(ls('P*'));
        length_choices = length(e.P.cur.choice);
    end
    
    for ii = 1 : length_choices
        disp(num2str(ii))
        if tworuns
            if ii <= length_run % check if still in the first run
                % first run
                % this is the current run
                nn                   = 1;
                time_corr_second_run = 0;
                % run length correction
                rlc                  = 0;
            else
                % second run
                nn                   = 2;
                time_corr_second_run = time_corr_second_run_fixed;
                % run length correction
                % subtract rlc from every ii such that we end up with two subcells of the same length for each run
                rlc = length_run;
                if strcmp(cur_sub,'VPPG0115')
                    cd(root_dir)
                    cd(cur_sub)
                    cd('Behav\PDT')
                    cd(p_structs{2})
                    load(ls('P*'))
                    time_corr_second_run = 0;
                end
            end
        else
            % single session (like a long first run)
            time_corr_second_run = 0;
            % run length correction
            rlc = 0;
        end
        
        % MISSING AS A CENSORING REGRESSOR
        trial{nn}{ii-rlc} = [ceil(((P.t.cur_trial.iti_on(ii)+time_corr_second_run)/tr)), floor(((P.t.cur_trial.stimLA_off(ii)+time_corr_second_run)/tr))];
        if ~(P.cur.choice(ii) > 0 && P.cur.choice(ii) < 5)
            % find limits of the trial -> don't cut out too much.
            % cut out images
            missing{nn}(trial{nn}{ii-rlc}(1):trial{nn}{ii-rlc}(2)) = 1;
            continue
        else
            missing{nn}(trial{nn}{ii-rlc}(1):trial{nn}{ii-rlc}(2)) = 0;
        end
        
        if P.cur.choice(ii) > 0 && P.cur.choice(ii) < 5
            % onsets for pic on
            onsets_cell{nn}(1)    = {[cell2mat(onsets_cell{nn}(1)) P.t.cur_trial.stim_on(ii) + ...
                time_corr_second_run]};
            durations_cell{nn}(1) = {[cell2mat(durations_cell{nn}(1)) P.t.cur_trial.stimLA_watch(ii) ...
                - P.t.cur_trial.stim_on(ii)]};
            
            for pm = 1:length(pic_pres)
                pmod_cell{nn}(1).param(pm) = {[pmod_cell{nn}(1).param{pm} 0]};
            end
            
            switch P.cur.cat(ii)
                
                case 1 % gambling pic
                    cp = 1; % current picture category
                    pmod_cell{nn}(1).param{cp}(length(pmod_cell{nn}(1).param{cp})) = 1;
                    
                case 2  % negative pic
                    cp = 2; % current picture category
                    pmod_cell{nn}(1).param{cp}(length(pmod_cell{nn}(1).param{cp})) = 1;
                    
                case 3  % positive pic
                    cp = 3; % current picture category
                    pmod_cell{nn}(1).param{cp}(length(pmod_cell{nn}(1).param{cp})) = 1;
                case 6
                    cp = 4;
                case 5
                    cp = 5; % missing; not doing anything
                otherwise
                    error('a cat case that is unaccounted for!')
            end
            
            % PIC PLUS GAMBLE ON onsets and durations
            onsets_cell{nn}(2)    = {[cell2mat(onsets_cell{nn}(2)) P.t.cur_trial.stimLA_watch(ii) ...
                + time_corr_second_run]};
            durations_cell{nn}(2) = {[cell2mat(durations_cell{nn}(2)) P.t.cur_trial.stimLA_on(ii) ...
                - P.t.cur_trial.stimLA_watch(ii)]};
            
            % gain
            cur_gain     = agk_recode(str2double(cell2mat(P.gain.strings(P.cur.gamble{ii}(1)))),osg,nsg) - mean_gain;
            gain_orig    = cur_gain + mean_gain;
            % loss (here changed loss to abs. loss)
            cur_loss     = agk_recode(abs(str2double(cell2mat(P.loss.strings(P.cur.gamble{ii}(2))))),osl,nsl) - abs(mean_loss);
            loss_orig    = cur_loss+abs(mean_loss);
            % ed
            cur_point    = [gain_orig;loss_orig;0];
            ed           = agk_get_ed(cur_point,sp,vec);
            
            % gain
            pmod_cell{nn}(2).param(1) = {[pmod_cell{nn}(2).param{1} cur_gain]};
            
            % loss (here changed loss to abs. loss)
            pmod_cell{nn}(2).param(2) = {[pmod_cell{nn}(2).param{2} cur_loss]};
            
            % ed
            pmod_cell{nn}(2).param(3) = {[pmod_cell{nn}(2).param{3} ed]};
            
            %keyboard
            if acc_rec
                % get the acceptance
                cur_acc = agk_recode(P.cur.choice(ii),[1,2,3,4,5],[1,1,0,0,5]);
                
                % acc
                pmod_cell{nn}(2).param(4)   = {[pmod_cell{nn}(2).param{4} cur_acc]}; % acc
                
                % accXpic
                pmod_cell{nn}(2).param(4+1) = {[pmod_cell{nn}(2).param{4+1} 0]};     % accXgam
                pmod_cell{nn}(2).param(4+2) = {[pmod_cell{nn}(2).param{4+2} 0]};     % accXneg
                pmod_cell{nn}(2).param(4+3) = {[pmod_cell{nn}(2).param{4+3} 0]};     % accXpos
                
                if cp < 4
                    pmod_cell{nn}(2).param{4+(cp)}(end) = cur_acc;                   % accXcur category
                end
                
%                 % accXpic
%                 pmod_cell{nn}(2).param(4)   = {[pmod_cell{nn}(2).param{4} 0]};      % gam
%                 pmod_cell{nn}(2).param(4+1) = {[pmod_cell{nn}(2).param{4+1} 0]};    % neg
%                 pmod_cell{nn}(2).param(4+2) = {[pmod_cell{nn}(2).param{4+2} 0]};    % pos
%                 pmod_cell{nn}(2).param(4+3) = {[pmod_cell{nn}(2).param{4+3} 0]};    % neu
%                 pmod_cell{nn}(2).param{4+(cp-1)}(end) = cur_acc;                    % cur category
            end
            
            % PIC PLUS GAMBLE ON PLUS OPTIONS ON
            onsets_cell{nn}(3)   = {[cell2mat(onsets_cell{nn}(3)) P.t.cur_trial.stimLA_on(ii) ...
                + time_corr_second_run]};
            durations_cell{nn}(3) = {[cell2mat(durations_cell{nn}(3)) P.t.cur_trial.stimLA_off(ii) ...
                - P.t.cur_trial.stimLA_on(ii)]};
            % NO PARAM MOD
            
        end
    end
    
    
    cd(fullfile(root_dir, cur_sub, 'MRT\NIFTI\PDT'))
    found_epi = 0;
    cur_epi_dirs  = cellstr(ls('*_epi*'));
    cur_MoCo_dirs = cellstr(ls('*_MoCo*'));
    if length(cur_epi_dirs{1}) > 0 % check if we're working with epis
        found_epi = 1;
        cur_fold = cur_epi_dirs;
    elseif length(cur_MoCo_dirs{1}) > 0 % check if we're working with MoCo
        found_epi = 1;
        cur_fold = cur_MoCo_dirs;       % store current folder name
    end
    found_epi = 1;
    
    cd(cur_fold{1})
    load(ls('Preprocessing*')); % for microtime resolution
    preproc_job = matlabbatch;
    cur_nslices   = preproc_job{1}.spm.temporal.st.nslices;  % store slices etc in batch
    cur_refslice  = preproc_job{1}.spm.temporal.st.refslice;
    
    % load template for ss model
    if acc_rec == 1
        load(cur_tmpl{2})
    elseif acc_rec == 0
        if tworuns
            load(cur_tmpl{3})
        else
            load(cur_tmpl{1})
        end
    end
    
    % Concatenate Missing and multiple regressor .txt file
    if physio_inc && exist(fullfile(root_dir, cur_sub, 'Physio'))
        % with physio correction available
        cd(fullfile(root_dir, cur_sub, 'Physio\results'))
        corr_m_r = cellstr(ls('corr*'));
        if isempty(corr_m_r{1})
            mult_reg = cellstr(ls('multiple_regressors*'));
            for mm = 1:num_run
                % padding missing with 0's to make up for post-experiment scans
                file_mult      = load(mult_reg{mm});
                to_pad         = length(file_mult(:,1)) - length(missing{mm});
                if to_pad < 0
                    missing{mm} = missing{mm}(1:end+to_pad);
                elseif to_pad > 0
                    missing{mm} = [missing{mm} zeros(1,to_pad)];
                end
                file_mult_corr = horzcat(file_mult, missing{mm}');
                save(['corr_multiple_regressors_run_',num2str(mm),'.txt'],'file_mult_corr', '-ascii')
            end
        end
    elseif physio_inc == 0 || ~exist(fullfile(root_dir, cur_sub, 'Physio'))
        % no physio correction
        for mm = 1:num_run
            disp(['No multiple regressor file for subject ', cur_sub,'. Missing trials are added to rp file'])
            cd(fullfile(root_dir, cur_sub, 'MRT\NIFTI\PDT',cur_fold{mm}))
            corr_rp = cellstr(ls('corr*'));
            if isempty(corr_rp{1})
                rp = load(ls('rp*'));
                % padding missing with 0's to make up for post-experiment scans
                to_pad = length(rp(:,1)) - length(missing{mm});
                if to_pad < 0 
                    missing{mm} = missing{mm}(1:end+to_pad);
                elseif to_pad > 0
                    missing{mm} = [missing{mm} zeros(1,to_pad)];
                end

                rp_corr = horzcat(rp, missing{mm}');
                save(['corr_rp_aepi_run',num2str(mm),'.txt'], 'rp_corr', '-ascii')
            end
        end
    end
    
    %% Fill the batch
    for nr = 1:num_run
        cd(root_dir);  % add path
        cd(cur_sub);
        cd('MRT/NIFTI/PDT')
        cd(cur_fold{nr})
        matlabbatch{1}.spm.stats.fmri_spec.sess(nr).scans = [];
        matlabbatch{1}.spm.stats.fmri_spec.sess(nr).scans = cellstr(spm_select('FPList',pwd,'swuaepi'));
        
        % if physio should be included, go into Physio results folder
        % and load .txt files from physio correction
        if physio_inc && exist(fullfile(root_dir, cur_sub, 'Physio'))
            cd(fullfile(root_dir, cur_sub, 'Physio\results'))
            % create multi_reg subfield
            mult_reg = cellstr(ls('corr*'));
            matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = [];
            matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = {fullfile(root_dir, cur_sub, 'Physio/results',mult_reg{nr})};
            
        else
            % else load rp file from preprocession
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = [];
            matlabbatch{1}.spm.stats.fmri_spec.sess(nr).multi_reg = cellstr(spm_select('FPList',pwd,'^corr'));
        end
        
        
        if found_epi == 0 % incase no scans were found: throw away
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {''};
        end
        
        cd(root_dir);  % add path
        cd(cur_sub);
        cd('MRT/NIFTI/PDT')
    end
    
    % create the results dir
    mkdir(pwd,'results')
    cd('results')
    
    if exist([pwd filesep results_ssdir])
        if ow == 1
            cmd_rmdir([pwd filesep results_ssdir])
        elseif ow == 0
            disp('Results dir already present. Overwrite not allowed, I will skip this subject.')
            error_message =  [results_ssdir ' ' cur_sub ' results dir already present. Skipped.'];
            return
        end
    end
    
    agk_mkdir_ex(pwd,results_ssdir)
    cd(results_ssdir)
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(pwd);
    
    if behav_data == 1
        for n = 1:num_run
            onsets = onsets_cell{n};
            durations = durations_cell{n};
            pmod = pmod_cell{n};
            
            if isequal(length(names), length(onsets), length(durations))
                save(['mult_cond_',num2str(n),'.mat'],'names','onsets','durations','pmod');
                matlabbatch{1}.spm.stats.fmri_spec.sess(n).multi = cellstr(...
                    spm_select('FPList',pwd,['mult_cond_',num2str(n),'.mat']));
            else
                disp(['names, onsets, durations for subject ', cur_sub, ' are of different size.'])
                error_msg = ['names, onsets, durations for subject ', cur_sub, ' are of different size.'];
                return
            end
            
        end
        
    end
    
    
    % fill in microtime resolution
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = cur_nslices;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = cur_refslice;
    % specifiy an explicit mask, which voxels shall only be analysed
    matlabbatch{1}.spm.stats.fmri_spec.mask = cellstr(expl_mask);
    % specify TR-time
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr;
    
    %% contrast manager
    % get the names
    all_con_names = {};
    ct            = 0;
    for kk = 1:length(names)
        ct                = ct + 1;
        all_con_names{ct} = names{kk};
        for ll = 1:length(pmod(kk).name)
            % params
            ct = ct + 1;
            all_con_names{ct} = [names{kk} 'x' pmod(kk).name{ll}];
        end
    end
    
    con.contrastNames   = all_con_names;
    con.contrastType    = cellstr(repmat('T',length(all_con_names),1));
    con.contrastWeights = agk_basic_t_cons(length(con.contrastType));
    % if tworuns
    %     con.contrastWeights = agk_basic_t_cons_2sess(length(con.contrastType),...
    %         length(con.contrastType));
    % end
    
    % replicate and scale
    con.contrastRep = cellstr(repmat('replsc',length(con.contrastNames),1))';
    
    for jj = 1:length(con.contrastType)
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.name    = con.contrastNames{jj};
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.convec  = con.contrastWeights{jj};
        matlabbatch{3}.spm.stats.con.consess{jj}.tcon.sessrep = con.contrastRep{jj};
    end
    
    if behav_data == 1
        save('design.mat','matlabbatch');
    end
    
    try
        if run_it == 1
            spm_jobman('run',matlabbatch);
        end
        error_message = [results_ssdir ' ' cur_sub ' Estimation successfull.'];
    catch MExc
        error_message = {MExc,[results_ssdir ' ' cur_sub ' Estimation not successfull.']};
        cd(root_dir)
        return
    end
    cd(root_dir);
catch MExc
    disp('Something went wrong');
    error_message = MExc;
    disp(MExc.message)
end

end