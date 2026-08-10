classdef ea_sweetspot < handle
    % Sweetspot analysis class to handle visualizations of sweetspots in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        M % content of lead group project
        resultfig % figure handle to plot results
        ID % name / ID of sweetspot object
        posvisible = 1 % sweetspot visible
        negvisible = 0 % sourspot visible

        efieldthreshold = 200
        statlevel = 'VTAs' % stats metric to use, 1 = active coordinates, 2 = efields, 3 = vtas
        stattest = 'T-Test';
        stat0hypothesis = 'Zero';
        statimpthreshold = 0;
        statNthreshold = 0;
        statamplitudecorrection = 'None';
        statnormalization = 'None';
        corrtype = 'Spearman' % correlation strategy in case of using E-Fields.
        coverthreshold = 20; % of vtas needed to cover a single voxel to be considered
        posBaseColor = [1,1,1] % positive main color
        posPeakColor = [0.9176,0.2000,0.1373] % positive peak color

        negBaseColor = [1,1,1] % negative main color
        negPeakColor = [0.2824,0.6157,0.9725] % negative peak color

        splitbygroup = 0
        showsignificantonly = 1
        alphalevel = 0.05
        multcompstrategy = 'Uncorrected'; % could be 'Bonferroni'

        autorefresh=1;

        results
        % Subfields:
        cvlivevisualize = 0; % if set to 1 shows crossvalidation results during processing.
        basepredictionon = 'Mean of Scores';
        spotdrawn % struct contains sweetspot drawn in the resultfig
        drawobject % actual streamtube handle
        patientselection % selected patients to include. Note that connected fibers are always sampled from all (& mirrored) VTAs of the lead group file
        setlabels={};
        setselections={};
        customselection % selected patients in the custom test list
        allpatients % list of all patients (as from M.patient.list)
        mirrorsides = 0 % flag to mirror VTAs / Efields to contralateral sides using ea_flip_lr_nonlinear()
        responsevar % response variable
        responsevarlabel % label of response variable
        covars = {} % covariates
        covarlabels = {} % covariate labels
        analysispath % where to store results
        leadgroup % redundancy protocol only, path to original lead group project
        colorbar % colorbar information
        % stats: (how many fibers available and shown etc for GUI)
        stats
        % additional settings:
        rngseed = 'default';
        Nperm = 1000 % how many permutations in leave-nothing-out permtest strategy
        permtestMaxWorkers = 3 % cap on parallel workers used by permtest/predpermtest's parfor loops. Each worker holds its own full copy of obj (incl. obj.results.efield), so too many workers can exhaust memory on large analyses -- tune per-machine if needed.
        predpermtestBatchSize = 250 % predpermtest processes permutations in batches of this size (reading from a slim, single-precision companion file via ea_sweetspot_permtest_slim) instead of loading all Nperm rows of Rperm/pperm at once -- lower this further on tighter-memory machines, at the cost of more, smaller disk reads.
        stratifyPermutationsByGroup = 0 % if true, permtest shuffles obj.responsevar only within obj.M.patient.group buckets. Independent of splitbygroup (which is mainly a visualization/coloring setting) -- does not affect or get affected by it.
        kfold = 5 % divide into k sets when doing k-fold CV
        Nsets = 5 % divide into N sets when doing Custom (random) set test
        adjustforgroups = 1 % adjust correlations for group effects
        ExternalModelFile = 'None'
        useExternalModel = false;
        visualizeExternalModel = 0;
    end

    properties (Access = private)
        switchedFromSpace=3 % if switching space, this will protocol where from
    end

    methods
        function obj=ea_sweetspot(analysispath) % class constructor
            if exist('analysispath', 'var') && ~isempty(analysispath)
                obj.analysispath = analysispath;
                [~, ID] = fileparts(obj.analysispath);
                obj.ID = ID;
            end
        end

        function initialize(obj,datapath,resultfig)
            D = load(datapath, '-mat');
            if isfield(D, 'M') % Lead Group analysis path loaded
                obj.M = D.M;
                obj.leadgroup = datapath;

                testID = obj.M.guid;
                ea_mkdir([fileparts(obj.leadgroup),filesep,'sweetspots',filesep]);
                id = 1;
                while exist([fileparts(obj.leadgroup),filesep,'sweetspots',filesep,testID,'.sweetspot'],'file')
                    testID = [obj.M.guid, '_', num2str(id)];
                    id = id + 1;
                end
                obj.ID = testID;
                obj.resultfig = resultfig;

                if isfield(obj.M,'pseudoM')
                    obj.allpatients = obj.M.ROI.list;
                    obj.patientselection = 1:length(obj.M.ROI.list);
                    obj.M.root = [fileparts(datapath),filesep];
                    %obj.M.patient.list=obj.M.ROI.list; % copies
                    obj.M.patient.list = cell(size(obj.M.ROI.list,1), 1);
                    for i = 1:size(obj.M.ROI.list,1)
                        obj.M.patient.list{i,1} = obj.M.ROI.list{i,1};
                    end
                    obj.M.patient.group=obj.M.ROI.group; % copies
                else
                    obj.allpatients = obj.M.patient.list;
                    obj.patientselection = obj.M.ui.listselect;
                end

                obj.responsevarlabel = obj.M.clinical.labels{1};
                obj.covarlabels={};
            elseif  isfield(D, 'sweetspot')  % Saved sweetspot class loaded
                props = properties(D.sweetspot);
                for p =  1:length(props) %copy all public properties
                    if ~(strcmp(props{p}, 'analysispath') && ~isempty(obj.analysispath) ...
                            || strcmp(props{p}, 'ID') && ~isempty(obj.ID))
                        obj.(props{p}) = D.tractset.(props{p});
                    end
                end
                clear D
            else
                ea_error('You have opened a file of unknown type.')
                return
            end
            obj.calculate;
        end

        function calculate(obj)
            % in case of the sweetspot explorer, calculate rather means to
            % gather all E-Fields. To keep consistency of the logic with
            % discfiberexplorer and networkmappingexplorer, we will keep
            % the same name (calculate) for the function, nonetheless.

            % check that results aren't already there
            if ~isempty(obj.results) % vtas already gathered in
                return
            end

            if isfield(obj.M,'pseudoM')
                vatlist = obj.M.ROI.list;
            else
                vatlist = ea_sweetspot_getvats(obj);
            end
            [AllX,space] = ea_exportefieldmapping(vatlist,obj);

            obj.results.efield = AllX;
            obj.results.space = space;

            if ~isfield(obj.M,'pseudoM')
                % get active coordinates, as well
                for pt=1:length(obj.M.patient.list)
                    for side=1:2
                        obj.results.activecnt{side}(pt,:)=...
                            mean(obj.M.elstruct(pt).coords_mm{side}(find(obj.M.S(pt).activecontacts{side}),:),1); %#ok<FNDSB> % find is necessary here
                    end
                end
                for side=1:2
                    obj.results.activecnt{side}=[obj.results.activecnt{side};ea_flip_lr_nonlinear(obj.results.activecnt{side})];
                end
            end

        end

        function Amps = getstimamp(obj)
            Amps=zeros(length(obj.M.patient.list),2);
            for pt=1:length(obj.M.patient.list)
                for side=1:2
                    thisamp=obj.M.stats(pt).ea_stats.stimulation.vat(side).amp;
                    thisamp(thisamp==0)=nan;
                    Amps(pt,side)=ea_nanmean(thisamp');
                end
            end
        end

        function VTAvolumes = getvtavolumes(obj)
            if ~isfield(obj.M.stats(1).ea_stats.stimulation.vat(1),'volume')
               VTAvolumes = obj.getstimamp;
                warning('No VTA volumes found. Using stimulation amplitudes instead. Re-run stats in Lead-group to obtain volumes.');
                return
            end
            VTAvolumes=zeros(length(obj.M.patient.list),2);
            for pt=1:length(obj.M.patient.list)
                for side=1:2
                    VTAvolumes(pt,side)=obj.M.stats(pt).ea_stats.stimulation.vat(side).volume;
                end
            end
        end

        function Efieldmags = getefieldmagnitudes(obj)
            if ~isfield(obj.M.stats(1).ea_stats.stimulation.efield(1),'volume')
                Efieldmags = obj.getstimamp;
                warning('No Efield magnitude sums found. Using stimulation amplitudes instead. Re-run stats in Lead-group to obtain values.');
                return
            end
            Efieldmags=zeros(length(obj.M.patient.list),2);
            for pt=1:length(obj.M.patient.list)
                for side=1:2
                    try
                        if isempty(obj.M.stats(pt).ea_stats.stimulation.efield(side).volume)
                            val=0;
                        else
                            val=obj.M.stats(pt).ea_stats.stimulation.efield(side).volume;
                        end
                        Efieldmags(pt,side)=val;
                    catch % could be efield(side) is not defined.
                        Efieldmags(pt,side)=0;
                    end
                end
            end
        end

        function refreshlg(obj)
            if ~exist(obj.leadgroup,'file')
               msgbox('Groupan alysis file has vanished. Please select file.');
               [fn,pth]=uigetfile();
               obj.leadgroup=fullfile(pth,fn);
            end
            D = load(obj.leadgroup);
            obj.M = D.M;
        end

        function coh = getcohortregressor(obj)
            coh=ea_cohortregressor(obj.M.patient.group(obj.patientselection));
        end

        function [I, Ihat] = loocv(obj)
            rng(obj.rngseed);
            cvp = cvpartition(length(obj.patientselection), 'LeaveOut');
            [I, Ihat] = crossval(obj, cvp);
        end

        function [I, Ihat] = lococv(obj)
            if length(unique(obj.M.patient.group(obj.patientselection))) == 1
                ea_error(sprintf(['Only one cohort in the analysis.\n', ...
                    'Leave-One-Cohort-Out cross-validation not possible.']));
            end
            [I, Ihat] = crossval(obj, obj.M.patient.group(obj.patientselection));
        end

        function [I, Ihat] = kfoldcv(obj)
            rng(obj.rngseed);
            cvp = cvpartition(length(obj.patientselection), 'KFold', obj.kfold);
            [I, Ihat] = crossval(obj, cvp);
        end

        function [I, Ihat] = lno(obj, Iperm)
            rng(obj.rngseed);
            cvp = cvpartition(length(obj.patientselection), 'resubstitution');
            if ~exist('Iperm', 'var')
                [I, Ihat] = crossval(obj, cvp);
            else
                [I, Ihat] = crossval(obj, cvp, Iperm);
            end
        end

        function [I, Ihat] = crossval(obj, cvp, Iperm)
            if isnumeric(cvp) % cvp is crossvalind
                cvIndices = cvp;
                cvID = unique(cvIndices);
                cvp = struct;
                cvp.NumTestSets = length(cvID);
                for i=1:cvp.NumTestSets
                    cvp.training{i} = cvIndices~=cvID(i);
                    cvp.test{i} = cvIndices==cvID(i);
                end
            end

            % Check if patients are selected in the custom training/test list
            if isempty(obj.customselection)
                patientsel = obj.patientselection;
            else
                patientsel = obj.customselection;
            end

            if ~exist('Iperm', 'var') || isempty(Iperm)
                I = obj.responsevar(patientsel,:);
            else
                I = Iperm(patientsel,:);
            end

            % Ihat is the estimate of improvements (not scaled to real improvements)
            Ihat = nan(length(patientsel),2);

            for c=1:cvp.NumTestSets
                if cvp.NumTestSets ~= 1
                    fprintf(['\nIterating set: %0',num2str(numel(num2str(cvp.NumTestSets))),'d/%d\n'], c, cvp.NumTestSets);
                end

                if isobject(cvp)
                    training = cvp.training(c);
                    test = cvp.test(c);
                elseif isstruct(cvp)
                    training = cvp.training{c};
                    test = cvp.test{c};
                end

                if obj.useExternalModel == true && ~strcmp(obj.ExternalModelFile, 'None')
                    % load external model, and assign vals from the
                    % external model.
                    S=ea_sweetspot_importedModel2Efields(obj, obj.ExternalModelFile);;
                    if obj.cvlivevisualize
                        [vals] = S.model_vals;
                        obj.draw(vals);
                        drawnow;
                    else
                        [vals] = S.model_vals;

                    end
                else
                    if ~exist('Iperm', 'var')
                        if obj.cvlivevisualize
                            [vals] = ea_sweetspot_calcstats(obj, patientsel(training));
                            obj.draw(vals);
                            drawnow;
                        else
                            [vals] = ea_sweetspot_calcstats(obj, patientsel(training));
                        end
                    else
                        if obj.cvlivevisualize
                            [vals] = ea_sweetspot_calcstats(obj, patientsel(training), Iperm);
                            obj.draw(vals);
                            drawnow;
                        else
                            [vals] = ea_sweetspot_calcstats(obj, patientsel(training), Iperm);
                        end
                    end
                end
                for side=1:numel(vals)
                    if ~isempty(vals{1,side})
                        switch obj.statlevel % also differentiate between methods in the prediction part.
                            case 'VTAs'
                                efield = obj.results.efield{side}(patientsel(test),:)';
                                efield(~isnan(efield)) = efield(~isnan(efield)) > obj.efieldthreshold;
                                switch lower(obj.basepredictionon)
                                    case 'mean of scores'
                                        Ihat(test,side) = ea_nanmean(obj.maskvals(vals{1,side},obj.posvisible,obj.negvisible).*efield,1);
                                    case 'sum of scores'
                                        Ihat(test,side) = ea_nansum(obj.maskvals(vals{1,side},obj.posvisible,obj.negvisible).*efield,1);
                                    case 'peak of scores'
                                        Ihat(test,side) = ea_discfibers_getpeak(vals{1,side}.*efield, obj.posvisible, obj.negvisible, 'peak');
                                    case 'peak 5% of scores'
                                        Ihat(test,side) = ea_discfibers_getpeak(vals{1,side}.*efield, obj.posvisible, obj.negvisible, 'peak5');
                                end
                            case 'E-Fields'
                                switch lower(obj.basepredictionon)
                                    case 'profile of scores: spearman'
                                        Ihat(test,side) = atanh(ea_corr(obj.maskvals(vals{1,side},obj.posvisible,obj.negvisible),obj.results.efield{side}(patientsel(test),:)','spearman'));
                                    case 'profile of scores: pearson'
                                        Ihat(test,side) = atanh(ea_corr(obj.maskvals(vals{1,side},obj.posvisible,obj.negvisible),obj.results.efield{side}(patientsel(test),:)','pearson'));
                                   case 'profile of scores: bend'
                                        Ihat(test,side) = atanh(ea_corr(obj.maskvals(vals{1,side},obj.posvisible,obj.negvisible),obj.results.efield{side}(patientsel(test),:)','bend'));
                                    case 'mean of scores'
                                        Ihat(test,side) = ea_nanmean(obj.maskvals(vals{1,side},obj.posvisible,obj.negvisible).*obj.results.efield{side}(patientsel(test),:)',1);
                                    case 'sum of scores'
                                        Ihat(test,side) = ea_nansum(obj.maskvals(vals{1,side},obj.posvisible,obj.negvisible).*obj.results.efield{side}(patientsel(test),:)',1);
                                    case 'peak of scores'
                                        Ihat(test,side) = ea_discfibers_getpeak(vals{1,side}.*obj.results.efield{side}(patientsel(test),:)', obj.posvisible, obj.negvisible, 'peak');
                                    case 'peak 5% of scores'
                                        Ihat(test,side) = ea_discfibers_getpeak(vals{1,side}.*obj.results.efield{side}(patientsel(test),:)', obj.posvisible, obj.negvisible, 'peak5');
                                end
                        end
                    end
                end
            end

            % restore original view in case of live drawing
            if obj.cvlivevisualize
                obj.draw;
            end

            if cvp.NumTestSets == 1
                Ihat = Ihat(test,:);
                I = I(test);
            end

            if size(obj.responsevar,2)==2 % hemiscores
                Ihat = Ihat(:); % compare hemiscores (electrode wise)
                I = I(:);
            else
                Ihat = ea_nanmean(Ihat,2); % compare bodyscores (patient wise)
                % Ihat = Ihat(:,2); % test a single side
                % % option 1
                % for i = 1:length(Ihat)
                %     if abs(Ihat(i,1))>abs(Ihat(i,2))
                %         Ihat_new(i) = Ihat(i,1);
                %     else
                %         Ihat_new(i) = Ihat(i,2);
                %     end
                % end       
                %option 2
                % for i = 1:length(Ihat)
                %     if (Ihat(i,1) <0) && (Ihat(i,2) < 0)
                %         Ihat_new(i) = min(Ihat(i,:));
                %     else
                %         Ihat_new(i) = max(Ihat(i,:));
                %     end
                % end     
                % 
                % Ihat = Ihat_new';
            end
        end

        function [Iperm, Ihat, R0, R1, pperm, Rp95] = lnopb(obj, corrType)
            if ~exist('corrType', 'var')
                corrType = 'Spearman';
            end

            numPerm = obj.Nperm;

            Iperm = ea_shuffle(obj.responsevar, numPerm, obj.patientselection, obj.rngseed)';
            Iperm = [obj.responsevar, Iperm];
            Ihat = cell(numPerm+1, 1);

            R = zeros(numPerm+1, 1);

            for perm=1:numPerm+1
                if perm==1
                    fprintf('Calculating without permutation\n\n');
                    [~, Ihat{perm}] = lno(obj);
                else
                    fprintf('Calculating permutation: %d/%d\n\n', perm-1, numPerm);
                    [~, Ihat{perm}] = lno(obj, Iperm(:, perm));
                end

                R(perm) = corr(Iperm(obj.patientselection,perm),Ihat{perm},'type',corrType,'rows','pairwise');
            end

            % generate null distribution
            R1 = R(1);
            R0 = sort(abs(R(2:end)),'descend');
            Rp95 = R0(round(0.05*numPerm));
            pperm = mean(abs(R0)>=abs(R1));
            disp(['Permuted p = ',sprintf('%0.2f',pperm),'.']);

            % Return only selected I
            Iperm = Iperm(obj.patientselection,:);
        end

        function permresults = permtest(obj, Nperm, corrType)
            % Voxelwise permutation test of the sweetspot correlation map.
            % Permutes obj.responsevar (group-restricted, NaN-excluded, tied
            % across sides/mirrors since permutation happens before the L/R &
            % mirror expansion in ea_sweetspot_calcstats), rebuilds the map
            % for each permutation, and saves the resulting empirical +
            % permuted voxelwise R/p maps for later statistical comparison.
            if ~exist('Nperm', 'var') || isempty(Nperm)
                Nperm = obj.Nperm;
            end
            if ~exist('corrType', 'var') || isempty(corrType)
                corrType = obj.corrtype;
            end

            if ~strcmp(obj.statlevel, 'E-Fields') || ~strcmp(obj.stattest, 'Correlations')
                ea_error('permtest is currently only implemented for statlevel = ''E-Fields'' with stattest = ''Correlations''.');
            end

            if size(obj.responsevar, 2) > 1
                ea_error('Hemiscore responsevar (2 columns) is not yet supported for permutation testing.');
            end

            patsel = obj.patientselection;

            if obj.stratifyPermutationsByGroup
                groupvec = obj.M.patient.group;
            else
                groupvec = [];
            end

            [Iperm, PermIdx] = ea_shuffle_grouped(obj.responsevar, Nperm, patsel, groupvec, obj.rngseed);

            % Descriptive base filename, also used as this run's dedicated folder
            % name -- encodes the options that most affect the result, plus a
            % timestamp so repeated runs never silently collide/overwrite each
            % other. Every output this run produces (the saved .mat, nifti
            % exports, README) lives together in sweetspots/<basefname>/, not
            % flat alongside every other run's files.
            if obj.stratifyPermutationsByGroup
                stratTag = 'stratified';
            else
                stratTag = 'pooled';
            end
            if obj.mirrorsides
                mirrorTag = 'mirrored';
            else
                mirrorTag = 'nonmirrored';
            end
            basefname = sprintf('%s_permtest_N%d_%s_%s_%s_%s', obj.ID, Nperm, stratTag, mirrorTag, corrType, datestr(now, 'yyyymmdd_HHMMSS'));

            outdir = [fileparts(obj.leadgroup), filesep, 'sweetspots', filesep, basefname, filesep];
            suffix = 1;
            while exist(outdir, 'dir') % timestamp makes this vanishingly rare, but never silently merge into an existing run's folder
                suffix = suffix + 1;
                outdir = [fileparts(obj.leadgroup), filesep, 'sweetspots', filesep, basefname, '_', num2str(suffix), filesep];
            end
            ea_mkdir(outdir);

            fprintf('Calculating empirical (unpermuted) sweetspot map...\n');
            [Remp, pemp, gvalFixed, gpatselFixed, thisvalsFixed, nanidxFixed] = ea_sweetspot_calcstats(obj, patsel, obj.responsevar, true); % skipsigthresh=true: always store raw, unmasked values
            % gvalFixed/gpatselFixed (coverage-masked efield data & patient
            % selection) and thisvalsFixed/nanidxFixed (the patient-sliced,
            % NaN-filtered correlation input derived from them) don't depend on
            % I/Iperm -- reused for every permutation below instead of being
            % recomputed (which forces full-matrix copies) on every one of the
            % Nperm calls.

            % export the empirical (raw, unmasked) R-map to nifti by default, so it
            % can be sanity-checked against any map generated the usual way (e.g. obj.draw()).
            for gi = 1:size(Remp,1)
                if all(cellfun(@isempty, Remp(gi,:)))
                    continue
                end
                ea_sweetspot_vals2nii(obj.results.space, Remp(gi,:), outdir, sprintf('%s_Remp_group%d', basefname, gi));
            end

            obj.capParpool;

            % ngroups/nsides are already knowable from Remp -- computed here,
            % before the parfor loop, so the loop itself can write directly
            % into its final, flat, single-precision form for the common
            % (single-group) case, rather than collecting every permutation's
            % full double-precision result into an Rrow/prow cell array first.
            % That intermediate held EVERY permutation's full result in memory
            % AT ONCE, in double precision, the moment the parfor loop
            % finished -- there is no way to reduce that after the fact by
            % clearing pieces of it sooner, since it's already fully built by
            % then. At Nperm=5000 on a real analysis this was ~123GB by
            % itself, which is what actually pushed a real run to ~180GB even
            % after adding row-by-row clearing (a real, but insufficient,
            % previous fix). Writing directly into preallocated
            % single-precision sliced parfor outputs avoids that intermediate
            % ever existing at all.
            ngroups = size(Remp,1);
            nsides = size(Remp,2);

            progress = 0;
            dq = parallel.pool.DataQueue;
            afterEach(dq, @(~) reportProgress());

            fprintf('Running %d permutations...\n', Nperm);

            permresults.Remp = Remp;
            permresults.pemp = pemp;

            if ngroups == 1
                % Fast path: flat, single-precision, top-level per-side
                % variables -- supports true row-chunked matfile() reads
                % downstream (predpermtest, ea_sweetspot_permtest_stats),
                % unlike a Rperm{side} cell of double arrays, which
                % matfile() can only load in full per element. This is the
                % layout ea_sweetspot_permtest_slim used to have to produce
                % as a separate, crash-prone conversion step; writing it
                % directly here, from the parfor loop itself, skips that
                % conversion, and its crash risk, entirely. See the
                % short-circuit at the top of ea_sweetspot_permtest_slim.m.
                %
                % Rperm_side1/2 (etc.) are "sliced" parfor output variables:
                % preallocated before the loop, each worker writes only its
                % own row (Rperm_side1(p,:) = ...), and MATLAB ships back
                % just that row rather than requiring the whole array to be
                % reconstructed from a cell of per-permutation pieces
                % afterward.
                Rperm_side1 = nan(Nperm, numel(Remp{1,1}), 'single');
                pperm_side1 = nan(Nperm, numel(Remp{1,1}), 'single');
                if nsides == 2
                    Rperm_side2 = nan(Nperm, numel(Remp{1,2}), 'single');
                    pperm_side2 = nan(Nperm, numel(Remp{1,2}), 'single');
                    parfor p = 1:Nperm
                        [v, pv] = ea_sweetspot_calcstats(obj, patsel, Iperm(:,p), true, gvalFixed, gpatselFixed, thisvalsFixed, nanidxFixed);
                        Rperm_side1(p,:) = single(v{1}');
                        pperm_side1(p,:) = single(pv{1}');
                        Rperm_side2(p,:) = single(v{2}');
                        pperm_side2(p,:) = single(pv{2}');
                        send(dq, 1);
                    end
                    permresults.Rperm_side2 = Rperm_side2;
                    permresults.pperm_side2 = pperm_side2;
                    clear Rperm_side2 pperm_side2
                else
                    parfor p = 1:Nperm
                        [v, pv] = ea_sweetspot_calcstats(obj, patsel, Iperm(:,p), true, gvalFixed, gpatselFixed, thisvalsFixed, nanidxFixed);
                        Rperm_side1(p,:) = single(v{1}');
                        pperm_side1(p,:) = single(pv{1}');
                        send(dq, 1);
                    end
                end
                permresults.Rperm_side1 = Rperm_side1;
                permresults.pperm_side1 = pperm_side1;
                clear Rperm_side1 pperm_side1
            else
                % Multiple groups (obj.splitbygroup): the fast flat layout
                % above -- and everything downstream that reads it
                % (predpermtest, ea_sweetspot_permtest_slim) -- has never
                % supported more than group 1, so this rarer path keeps
                % collecting every permutation's full {group,side} result via
                % Rrow/prow first, then reshaping into the original nested
                % {group,side} double-cell layout, unchanged -- preallocating
                % per-group/side sliced outputs the same way as above isn't
                % worth the complexity for this already-slower legacy path.
                % Each Rrow{p}{gi}/prow{p}{gi} is still freed as soon as it's
                % copied out below, to at least avoid holding it twice over.
                Rrow = cell(1, Nperm);
                prow = cell(1, Nperm);
                parfor p = 1:Nperm
                    [v, pv] = ea_sweetspot_calcstats(obj, patsel, Iperm(:,p), true, gvalFixed, gpatselFixed, thisvalsFixed, nanidxFixed);
                    Rrow{p} = v;
                    prow{p} = pv;
                    send(dq, 1);
                end

                Rperm = cell(size(Remp));
                pperm = cell(size(Remp));
                for gi = 1:numel(Remp)
                    if isempty(Remp{gi})
                        continue
                    end
                    Rperm{gi} = nan(Nperm, numel(Remp{gi}));
                    pperm{gi} = nan(Nperm, numel(Remp{gi}));
                    for p = 1:Nperm
                        Rperm{gi}(p,:) = Rrow{p}{gi}';
                        pperm{gi}(p,:) = prow{p}{gi}';
                        Rrow{p}{gi} = [];
                        prow{p}{gi} = [];
                    end
                end
                permresults.Rperm = Rperm;
                permresults.pperm = pperm;
                clear Rrow prow
            end

            permresults.PermIdx = PermIdx;
            permresults.Nperm = Nperm;
            permresults.rngseed = obj.rngseed;
            permresults.corrtype = corrType;
            permresults.splitbygroup = obj.splitbygroup; % visualization/coloring setting -- not what governed the shuffle, kept only for reference
            permresults.stratifyPermutationsByGroup = obj.stratifyPermutationsByGroup; % this is what actually governed the shuffle
            permresults.mirrorsides = obj.mirrorsides;
            permresults.patientselection = patsel;
            permresults.statlevel = obj.statlevel;
            permresults.stattest = obj.stattest;
            permresults.coverthreshold = obj.coverthreshold;
            permresults.efieldthreshold = obj.efieldthreshold;
            permresults.space = obj.results.space; % training grid geometry, needed for out-of-sample reslicing later
            permresults.ID = obj.ID;
            permresults.leadgroup = obj.leadgroup;

            % outdir is a freshly-created, just-uniquified folder (see above) --
            % nothing else could already have a file here, so no filename
            % collision check is needed at this level (unlike outdir itself).
            outfile = [outdir, basefname, '.mat'];
            permresults.savedfile = outfile;
            save(outfile, '-struct', 'permresults', '-v7.3');
            fprintf('Saved permutation results to %s\n', outfile);

            readmeBody = sprintf([ ...
                'Voxelwise permutation test of the sweetspot correlation map for %s.\n\n', ...
                'Settings used for this run:\n', ...
                '  Nperm                      = %d\n', ...
                '  corrType                   = %s\n', ...
                '  stratifyPermutationsByGroup = %d\n', ...
                '  mirrorsides                = %d\n', ...
                '  statlevel                  = %s\n', ...
                '  stattest                   = %s\n', ...
                '  coverthreshold             = %g%%\n', ...
                '  efieldthreshold            = %g\n', ...
                '  rngseed                    = %s\n\n', ...
                'Files in this folder:\n', ...
                '  %s.mat  - Remp/pemp (empirical R/p per voxel) and Rperm/pperm\n', ...
                '                 (null distributions -- flat, single-precision, chunk-readable\n', ...
                '                 for the common single-group case; see ea_sweetspot_permtest_slim\n', ...
                '                 for the legacy/multi-group nested-cell format), plus PermIdx\n', ...
                '                 (the exact patient-shuffle used per permutation) and the settings\n', ...
                '                 block above.\n', ...
                '  *_Remp_group*.nii - the empirical (real, unpermuted) correlation map, UNMASKED\n', ...
                '                 (raw R everywhere, no significance thresholding applied yet --\n', ...
                '                 that happens in a separate post-hoc step, see\n', ...
                '                 ea_sweetspot_permtest_stats, which adds its own files/section to\n', ...
                '                 this same folder/README once run).\n'], ...
                obj.ID, Nperm, corrType, obj.stratifyPermutationsByGroup, obj.mirrorsides, obj.statlevel, ...
                obj.stattest, obj.coverthreshold, obj.efieldthreshold, obj.rngseed, basefname);
            ea_sweetspot_readme_append(outdir, 'Voxelwise permutation test (permtest)', readmeBody);

            function reportProgress()
                % Called on the client (not inside the workers) each time a
                % worker finishes one permutation, via the DataQueue above.
                progress = progress + 1;
                step = max(1, round(Nperm/20)); % ~5% increments
                if mod(progress, step) == 0 || progress == Nperm
                    fprintf('Permutation progress: %d/%d (%.0f%%)\n', progress, Nperm, 100*progress/Nperm);
                end
            end
        end

        function save(obj)
            sweetspot=obj;
            pth = fileparts(sweetspot.leadgroup);
            sweetspot.analysispath=[pth,filesep,'sweetspots',filesep,obj.ID,'.sweetspot'];
            ea_mkdir([pth,filesep,'sweetspots']);
            rf=obj.resultfig; % need to stash fig handle for saving.
            rd=obj.drawobject; % need to stash handle of drawing before saving.
            try % could be figure is already closed.
                setappdata(rf,['dt_',sweetspot.ID],rd); % store handle of tract to figure.
            end
            sweetspot.resultfig=[]; % rm figure handle before saving.
            sweetspot.drawobject=[]; % rm drawobject.
            save(sweetspot.analysispath,'sweetspot','-v7.3');
            obj.resultfig=rf;
            obj.drawobject=rd;
        end

        function export=draw(obj,vals)
            if obj.useExternalModel == true && ~strcmp(obj.ExternalModelFile, 'None') && obj.visualizeExternalModel == 1
                % load external model, and assign vals from the
                % external model.
                S=ea_sweetspot_importedModel2Efields(obj, obj.ExternalModelFile);
                [vals] = S.model_vals;           
            elseif ~exist('vals','var')
                [vals]=ea_sweetspot_calcstats(obj);
            end
            obj.spotdrawn.vals=vals;

            obj.stats.pos.shown(1)=sum(vals{1,1}>0);
            obj.stats.neg.shown(1)=sum(vals{1,1}<0);

            set(0,'CurrentFigure',obj.resultfig);

            dogroups=size(vals,1)>1; % if color by groups is set will be positive.
            if ~isfield(obj.M,'groups')
                obj.M.groups.group=ones(length(obj.M.patient.list),1);
                obj.M.groups.color=ea_color_wes('all');
            end
            linecols=obj.M.groups.color;
            if isempty(obj.drawobject) % check if prior object has been stored
                obj.drawobject=getappdata(obj.resultfig,['dt_',obj.ID]); % store handle of tract to figure.
            end
            for s=1:numel(obj.drawobject)
                for ins=1:numel(obj.drawobject{s})
                    try delete(obj.drawobject{s}{ins}.toggleH); end
                    try delete(obj.drawobject{s}{ins}.patchH); end
                    try delete(obj.drawobject{s}{ins}); end
                end
            end
            obj.drawobject={};

            % reset colorbar
            obj.colorbar=[];
            if ~any([obj.posvisible,obj.negvisible])
                export=nan;
            end

            for group=1:size(vals,1) % vals will have 1x2 in case of bipolar drawing and Nx2 in case of group-based drawings (where only positives are shown).
                % Vertcat all values for colorbar construction
                allvals = vertcat(vals{group,:});
                if isempty(allvals) || all(isnan(allvals))
                    ea_cprintf('CmdWinWarnings', 'Empty or all-nan value found!\n');
                    continue;
                else
                    allvals(isnan(allvals)) = 0;
                end

                if obj.posvisible && all(allvals<=0)
                    obj.posvisible = 0;
                    fprintf('\n')
                    warning('off', 'backtrace');
                    warning('No positive values found, posvisible is set to 0 now!');
                    warning('on', 'backtrace');
                    fprintf('\n')
                end

                if obj.negvisible && all(allvals>=0)
                    obj.negvisible = 0;
                    fprintf('\n')
                    warning('off', 'backtrace');
                    warning('No negative values found, negvisible is set to 0 now!');
                    warning('on', 'backtrace');
                    fprintf('\n')
                end

                colormap(gray);
                gradientLevel = length(gray);

                if dogroups
                    if obj.posvisible && ~obj.negvisible
                        voxcmap{group} = ea_colorgradient(gradientLevel, obj.posBaseColor, linecols(group,:));
                    elseif ~obj.posvisible && obj.negvisible
                        voxcmap{group} = ea_colorgradient(gradientLevel, linecols(group,:), obj.negBaseColor);
                    else
                        warndlg(sprintf(['Please choose either "Show Positive Regions" or "Show Negative Regions".',...
                            '\nShow both positive and negative regions is not supported when "Color by Group Variable" is on.']));
                        return;
                    end
                else
                    if obj.posvisible && obj.negvisible
                        cmapLeft = ea_colorgradient(gradientLevel/2, obj.negPeakColor, obj.negBaseColor);
                        cmapRight = ea_colorgradient(gradientLevel/2, obj.posBaseColor, obj.posPeakColor);
                        voxcmap{group} = [cmapLeft;cmapRight];
                    elseif obj.posvisible
                        voxcmap{group} = ea_colorgradient(gradientLevel, obj.posBaseColor, obj.posPeakColor);
                    elseif obj.negvisible
                        voxcmap{group} = ea_colorgradient(gradientLevel, obj.negPeakColor, obj.negBaseColor);
                    end
                end

                for side=1:size(vals,2)
                    res=obj.results.space{side};
                    res.dt(1) = 16;
                    res.img(:)=nan;
                    % Plot voxels if any survived
                    if obj.posvisible
                        % plot positives:
                        posvox=res;
                        posvox.img(:)=0;
                        posvox.img(vals{group,side}>0)=vals{group,side}(vals{group,side}>0);

                        pobj.nii=posvox;
                        pobj.name='Positive';
                        pobj.niftiFilename='Positive.nii';
                        pobj.binary=0;
                        pobj.usesolidcolor=0;
                        pobj.color=obj.posPeakColor;
                        pobj.colormap=ea_colorgradient(gradientLevel, obj.posBaseColor, obj.posPeakColor);
                        pobj.smooth=10;
                        pobj.hullsimplify=0.5;
                        pobj.threshold=0;
                        obj.drawobject{group,side}{1}=ea_roi('Positive.nii',pobj);

                        res=posvox; % keep copy for export
                    end

                    if obj.negvisible
                        % plot negatives:
                        negvox=res;
                        negvox.img(:)=0;
                        negvox.img(vals{group,side}<0)=-vals{group,side}(vals{group,side}<0);

                        pobj.nii=negvox;
                        pobj.name='Negative';
                        pobj.niftiFilename='Negative.nii';
                        pobj.binary=0;
                        pobj.usesolidcolor=0;
                        pobj.color=obj.negPeakColor;
                        pobj.colormap=ea_colorgradient(gradientLevel, obj.negPeakColor, obj.negBaseColor);
                        pobj.smooth=10;
                        pobj.hullsimplify=0.5;
                        pobj.threshold=0;
                        obj.drawobject{group,side}{2}=ea_roi('Negative.nii',pobj);

                        res.img(:)=nansum([res.img(:),-negvox.img(:)],2); % keep copy for export.
                    end
                    res.img(res.img==0)=nan;
                    export{side}=res;
                end

                % Set colorbar tick positions and labels
                if ~isempty(allvals)
                    if obj.posvisible && obj.negvisible
                        tick{group} = [1, gradientLevel/2-10, gradientLevel/2+11, length(voxcmap{group})];
                        poscbvals = sort(allvals(allvals>0));
                        negcbvals = sort(allvals(allvals<0));
                        ticklabel{group} = [negcbvals(1), negcbvals(end), poscbvals(1), poscbvals(end)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    elseif obj.posvisible
                        tick{group} = [1, length(voxcmap{group})];
                        posvals = sort(allvals(allvals>0));
                        ticklabel{group} = [posvals(1), posvals(end)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    elseif obj.negvisible
                        tick{group} = [1, length(voxcmap{group})];
                        negvals = sort(allvals(allvals<0));
                        ticklabel{group} = [negvals(1), negvals(end)];
                        ticklabel{group} = arrayfun(@(x) num2str(x,'%.2f'), ticklabel{group}, 'Uni', 0);
                    end
                end
            end

            if ~exist('export','var') % all empty
                for side=1:size(vals,2)
                    res=obj.results.space{side};
                    res.dt(1) = 16;
                    res.img(:)=nan;
                    export{side}=res;
                end
            end

            setappdata(obj.resultfig,['dt_',obj.ID],obj.drawobject); % store handle of surf to figure.

            % store colorbar in object
            if exist('voxcmap','var')
                setappdata(obj.resultfig, ['voxcmap',obj.ID], voxcmap);
                obj.colorbar.cmap = voxcmap;
                obj.colorbar.tick = tick;
                obj.colorbar.ticklabel = ticklabel;
            end
        end

        function predresults = predpermtest(obj, trainPermtestFile, sigMode, tail)
            % Out-of-sample validation of a voxelwise permutation test (built
            % via permtest() on a *different*, training ea_sweetspot object)
            % against this (test) object's own, never-permuted patient
            % scores. Each permuted (and the empirical) training map is
            % reprojected onto this object's own efield grid via a one-time
            % nearest-neighbor index map -- the geometric correspondence
            % only depends on the two objects' obj.results.space, not on the
            % map values, so it is computed once and reused for every
            % permutation instead of reslicing Nperm+1 times.
            %
            % sigMode controls uncorrected significance thresholding (using
            % obj.alphalevel, evaluated on this -- the test -- object) before
            % prediction, matching how the empirical model is normally
            % restricted to significant voxels only:
            %   'None'        - no thresholding, raw R everywhere (default)
            %   'Independent' - each map (empirical & every permutation) is
            %                   thresholded using its own p-values
            %   'Fixed'       - the empirical map's significant-voxel mask is
            %                   computed once and applied to every permutation
            %
            % tail controls which direction(s) of the null distribution count
            % as "as extreme as" the empirical prediction, for exceedCount/pperm:
            %   'right' (default) - Rpredperm >= Rpredemp (a hypothesis that
            %                       the map should only ever predict IMPROVEMENT)
            %   'left'             - Rpredperm <= Rpredemp
            %   'both'             - abs(Rpredperm) >= abs(Rpredemp); this
            %                        function's only behavior before tail existed

            if size(obj.responsevar, 2) > 1
                ea_error('Hemiscore responsevar (2 columns) is not yet supported for permutation testing.');
            end

            if ~exist('sigMode', 'var') || isempty(sigMode)
                sigMode = 'None';
            end
            if ~ismember(sigMode, {'None', 'Independent', 'Fixed'})
                ea_error('sigMode must be ''None'', ''Independent'', or ''Fixed''.');
            end

            if ~exist('tail', 'var') || isempty(tail)
                tail = 'right';
            end
            if ~ismember(tail, {'right', 'left', 'both'})
                ea_error('tail must be ''right'', ''left'', or ''both''.');
            end

            % Rperm/pperm are the memory-heavy part of a saved permtest file
            % (Nperm x Nvoxels double -- e.g. ~60GB combined across sides at
            % Nperm=5000 on a ~769K-voxel bilateral grid). Loading them in one
            % shot via a blanket load() has crashed MATLAB outright on a
            % 48GB-RAM machine, even in a fresh session with nothing else
            % running. ea_sweetspot_permtest_slim converts them (once --
            % cached on disk afterwards, so repeat calls against the same
            % training file skip straight past this) into flat,
            % single-precision, top-level per-side variables that support
            % true row-chunked reads via matfile, processed below in batches
            % of obj.predpermtestBatchSize permutations so peak memory stays
            % a small fraction of the full Nperm x Nvoxels size.
            slimFile = ea_sweetspot_permtest_slim(trainPermtestFile);
            trainMeta = load(slimFile, 'space', 'Remp', 'pemp', 'Nperm', 'corrtype', 'PermIdx', 'ID');
            mSlim = matfile(slimFile);

            nsides = numel(obj.results.space);
            srcIdx = cell(1, nsides);
            testUncovered = zeros(1, nsides);
            testTotal = zeros(1, nsides);
            trainUncovered = zeros(1, nsides);
            trainTotal = zeros(1, nsides);

            fprintf('Voxel coverage report (training grid vs. this object''s test grid):\n');
            for side = 1:nsides
                srcIdx{side} = ea_sweetspot_nnindexmap(trainMeta.space{side}, obj.results.space{side});

                testTotal(side) = numel(srcIdx{side});
                testUncovered(side) = sum(isnan(srcIdx{side}));

                trainTotal(side) = numel(trainMeta.space{side}.img);
                covered = unique(srcIdx{side}(~isnan(srcIdx{side})));
                trainUncovered(side) = trainTotal(side) - numel(covered);

                fprintf(['  Side %d: %d/%d (%.1f%%) test-grid voxels have no corresponding training voxel.\n', ...
                    '           %d/%d (%.1f%%) training-grid voxels are not represented anywhere in the test grid.\n'], ...
                    side, testUncovered(side), testTotal(side), 100*testUncovered(side)/testTotal(side), ...
                    trainUncovered(side), trainTotal(side), 100*trainUncovered(side)/trainTotal(side));
            end

            if isempty(obj.customselection)
                patsel = obj.patientselection;
            else
                patsel = obj.customselection;
            end
            Nperm = trainMeta.Nperm;

            % Uncorrected significance thresholding (training-grid space, before
            % reprojection -- Remp/pemp and Rperm/pperm share the same voxel
            % indexing, so this is a plain elementwise mask).
            alphalevel = obj.alphalevel;
            Remp_train = trainMeta.Remp(1,:);
            if ~strcmp(sigMode, 'None')
                for side = 1:nsides
                    nonsig = isnan(trainMeta.pemp{1,side}) | trainMeta.pemp{1,side} > alphalevel;
                    Remp_train{side}(nonsig) = nan;
                end
            end
            fixedMask = cell(1, nsides); % only used when sigMode == 'Fixed'
            if strcmp(sigMode, 'Fixed')
                for side = 1:nsides
                    fixedMask{side} = ~isnan(Remp_train{side}); % logical, training-grid space
                end
            end

            % validMask/efieldT are both independent of which map (empirical or
            % which permutation) is being predicted -- precomputed once and reused
            % below instead of being recomputed on every one of the Nperm+1 calls.
            validMask = cell(1, nsides);
            efieldT = cell(1, nsides);
            for side = 1:nsides
                validMask{side} = ~isnan(srcIdx{side});
                efieldT{side} = obj.results.efield{side}(patsel,:)';
            end

            fprintf('Predicting from empirical (unpermuted) training map...\n');
            empvals = cell(1, nsides);
            for side = 1:nsides
                empvals{side} = nan(numel(srcIdx{side}), 1);
                empvals{side}(validMask{side}) = Remp_train{side}(srcIdx{side}(validMask{side}));
            end
            [Rpredemp, Rpredemp_pval, Ihatemp] = obj.predictfromvals(empvals, patsel, trainMeta.corrtype, efieldT);
            if isnan(Rpredemp)
                ea_error('Empirical out-of-sample prediction is NaN -- cannot rank against the null distribution. Check basepredictionon/posvisible/negvisible settings and voxel coverage.');
            end
            fprintf('Empirical out-of-sample prediction: Rpredemp = %.4f (parametric p = %.4g).\n', Rpredemp, Rpredemp_pval);
            fprintf('Check this against your independently-computed empirical R now -- if it does not match, stop here (Ctrl+C) before the %d-permutation null distribution runs.\n', Nperm);

            batchSize = obj.predpermtestBatchSize;
            fprintf('Predicting from %d permuted training maps, in batches of %d (serial -- no parfor)...\n', Nperm, batchSize);
            Rpredperm = nan(Nperm, 1);
            corrType = trainMeta.corrtype;

            % Serial, not parfor: a persistent worker pool driven by dozens of
            % sequential parfor calls (one per batch) has been observed to
            % accumulate memory across calls until a worker gets OOM-killed
            % and the whole pool fails to recover mid-run -- a different,
            % harder-to-fix failure mode than the original blanket-load()
            % crash the batching itself solves. Running serially removes the
            % worker pool from the picture entirely, at the cost of using one
            % core instead of obj.permtestMaxWorkers.
            for batchStart = 1:batchSize:Nperm
                batchIdx = batchStart:min(batchStart+batchSize-1, Nperm);
                nb = numel(batchIdx);

                % Only this batch's rows ever touch memory -- mSlim.Rperm_sideN
                % is a plain top-level array (not nested in a cell), so matfile
                % reads exactly these rows from disk instead of materializing
                % the full Nperm x Nvoxels array.
                Rbatch1 = mSlim.Rperm_side1(batchIdx,:);
                pbatch1 = mSlim.pperm_side1(batchIdx,:);
                if nsides >= 2
                    Rbatch2 = mSlim.Rperm_side2(batchIdx,:);
                    pbatch2 = mSlim.pperm_side2(batchIdx,:);
                else
                    Rbatch2 = [];
                    pbatch2 = [];
                end

                RpredBatch = nan(nb, 1);
                for bi = 1:nb
                    permvals = cell(1, nsides);
                    for side = 1:nsides
                        permvals{side} = nan(numel(srcIdx{side}), 1);
                        valid = validMask{side};
                        if side == 1
                            row = double(Rbatch1(bi,:));
                            prow = double(pbatch1(bi,:));
                        else
                            row = double(Rbatch2(bi,:));
                            prow = double(pbatch2(bi,:));
                        end
                        switch sigMode
                            case 'Independent'
                                row(isnan(prow) | prow > alphalevel) = nan;
                            case 'Fixed'
                                row(~fixedMask{side}) = nan;
                        end
                        permvals{side}(valid) = row(srcIdx{side}(valid));
                    end
                    RpredBatch(bi) = obj.predictfromvals(permvals, patsel, corrType, efieldT);
                end
                Rpredperm(batchIdx) = RpredBatch;

                fprintf('  Completed permutations %d-%d of %d\n', batchIdx(1), batchIdx(end), Nperm);
            end

            % A permutation whose map produced no defined prediction (NaN -- e.g.
            % zero significant voxels under sigMode='Independent') is itself
            % evidence against the null being able to predict, not missing data:
            % it stays in the denominator and counts as not exceeding the
            % empirical result (comparisons against NaN are always false in
            % MATLAB, which already gives this behavior -- made explicit here
            % rather than left implicit, and reported instead of silent).
            nNaNperm = sum(isnan(Rpredperm));
            switch tail
                case 'right'
                    exceedCount = sum(Rpredperm >= Rpredemp);
                case 'left'
                    exceedCount = sum(Rpredperm <= Rpredemp);
                case 'both'
                    exceedCount = sum(abs(Rpredperm) >= abs(Rpredemp));
            end
            Rp0 = sort(abs(Rpredperm), 'descend'); % NaNs sort to the end automatically; unaffected by tail, always the two-tailed 95th-percentile diagnostic
            Rp95 = Rp0(round(0.05*Nperm));
            pperm = exceedCount / Nperm;
            fprintf('%d/%d permutations (%.1f%%) produced no defined prediction (NaN) -- counted as not exceeding the empirical result.\n', nNaNperm, Nperm, 100*nNaNperm/Nperm);
            disp(['Out-of-sample permuted p (tail=''', tail, ''') = ', sprintf('%0.3f', pperm), ' (empirical R ranks ', num2str(exceedCount), ' of ', num2str(Nperm), ').']);

            predresults.Rpredemp = Rpredemp;
            predresults.Rpredemp_pval = Rpredemp_pval; % parametric p-value of the empirical correlation itself (distinct from pperm, the permutation-based null p-value)
            predresults.Ihatemp = Ihatemp; % per-patient predicted score underlying Rpredemp
            predresults.Iemp = obj.responsevar(patsel); % this object's real, unpermuted scores -- both saved so the correlation plot can be regenerated standalone later (ea_sweetspot_predpermtest_plot), without a live object
            predresults.responsevarlabel = obj.responsevarlabel;
            predresults.Rpredperm = Rpredperm;
            predresults.nNaNperm = nNaNperm; % permutations with no defined prediction, counted as not exceeding Rpredemp (see pperm)
            predresults.exceedCount = exceedCount; % how many (of Nperm) permutations were as extreme as Rpredemp per `tail`; pperm = exceedCount/Nperm
            predresults.pperm = pperm;
            predresults.tail = tail;
            predresults.Rp95 = Rp95;
            predresults.trainPermIdx = trainMeta.PermIdx; % training cohort's shuffle indices (NOT the test cohort -- this object's own patients are never permuted). Columns correspond 1:1 to Rpredperm entries, traceable to the training permtest that produced each permuted map.
            predresults.trainPermtestFile = trainPermtestFile;
            predresults.trainID = trainMeta.ID;
            predresults.testID = obj.ID;
            predresults.corrtype = trainMeta.corrtype;
            predresults.basepredictionon = obj.basepredictionon;
            predresults.sigMode = sigMode;
            predresults.alphalevel = alphalevel;
            predresults.testUncoveredVoxels = testUncovered;
            predresults.testTotalVoxels = testTotal;
            predresults.trainUncoveredVoxels = trainUncovered;
            predresults.trainTotalVoxels = trainTotal;

            % Basefname (also used as this run's dedicated folder name) identifies
            % which training run this test cohort was checked against (a test
            % cohort may be validated against several different training analyses
            % for different purposes -- each gets its own folder instead of
            % overwriting the last one), plus sigMode and a timestamp.
            [~, trainBaseName] = fileparts(trainPermtestFile);
            basefname = sprintf('%s_predpermtest_vs_%s_%s_%s', obj.ID, trainBaseName, sigMode, datestr(now, 'yyyymmdd_HHMMSS'));

            outdir = [fileparts(obj.leadgroup), filesep, 'sweetspots', filesep, basefname, filesep];
            suffix = 1;
            while exist(outdir, 'dir') % timestamp makes this vanishingly rare, but never silently merge into an existing run's folder
                suffix = suffix + 1;
                outdir = [fileparts(obj.leadgroup), filesep, 'sweetspots', filesep, basefname, '_', num2str(suffix), filesep];
            end
            ea_mkdir(outdir);

            % outdir is a freshly-created, just-uniquified folder -- nothing else
            % could already have a file here, so no filename collision check is
            % needed at this level (unlike outdir itself, above).
            outfile = [outdir, basefname, '.mat'];
            predresults.savedfile = outfile;
            save(outfile, '-struct', 'predresults', '-v7.3');
            fprintf('Saved out-of-sample permutation prediction results to %s\n', outfile);

            readmeBody = sprintf([ ...
                'Out-of-sample validation of training permtest file:\n  %s\nagainst this (test) object''s (%s) own, never-permuted patient scores.\n\n', ...
                'Settings used for this run:\n', ...
                '  sigMode  = %s\n', ...
                '  tail     = %s\n', ...
                '  alphalevel = %g\n', ...
                '  basepredictionon = %s\n\n', ...
                'Result: Rpredemp = %.4f (parametric p = %.4g); out-of-sample permuted p = %.3f\n', ...
                '(empirical R ranks %d of %d permutations).\n\n', ...
                'Files in this folder:\n', ...
                '  %s.mat  - Rpredemp/Rpredperm (empirical + null out-of-sample predictions),\n', ...
                '                 Ihatemp/Iemp (per-patient predicted vs. real scores), and the\n', ...
                '                 settings/result summary above.\n'], ...
                trainPermtestFile, obj.ID, sigMode, tail, alphalevel, obj.basepredictionon, ...
                Rpredemp, Rpredemp_pval, pperm, exceedCount, Nperm, basefname);
            ea_sweetspot_readme_append(outdir, 'Out-of-sample prediction (predpermtest)', readmeBody);

            ea_sweetspot_predpermtest_plot(outfile); % shared with standalone re-plotting from a saved file
        end
    end

    methods (Access = private)
        function [Rpred, Rpred_p, Ihat] = predictfromvals(obj, vals, patsel, corrType, efieldTIn)
            % Replicates the E-Fields prediction logic in crossval() (see
            % obj.basepredictionon), but takes the voxelwise map(s) directly
            % rather than recomputing them, so it can be reused across many
            % permuted maps.
            %
            % efieldTIn (optional): precomputed obj.results.efield{side}(patsel,:)'
            % per side, since that transpose/slice is independent of vals and
            % otherwise gets recomputed on every one of the Nperm calls a caller
            % like predpermtest makes.
            if ~exist('efieldTIn', 'var')
                efieldTIn = {};
            end
            Ihat = nan(length(patsel), numel(vals));
            for side = 1:numel(vals)
                v = vals{side};
                if numel(efieldTIn) >= side && ~isempty(efieldTIn{side})
                    efieldT = efieldTIn{side};
                else
                    efieldT = obj.results.efield{side}(patsel,:)';
                end
                switch lower(obj.basepredictionon)
                    case 'profile of scores: spearman'
                        Ihat(:,side) = atanh(ea_corr(obj.maskvals(v,obj.posvisible,obj.negvisible), efieldT, 'spearman'));
                    case 'profile of scores: pearson'
                        Ihat(:,side) = atanh(ea_corr(obj.maskvals(v,obj.posvisible,obj.negvisible), efieldT, 'pearson'));
                    case 'profile of scores: bend'
                        Ihat(:,side) = atanh(ea_corr(obj.maskvals(v,obj.posvisible,obj.negvisible), efieldT, 'bend'));
                    case 'mean of scores'
                        Ihat(:,side) = ea_nanmean(obj.maskvals(v,obj.posvisible,obj.negvisible).*efieldT,1);
                    case 'sum of scores'
                        Ihat(:,side) = ea_nansum(obj.maskvals(v,obj.posvisible,obj.negvisible).*efieldT,1);
                    case 'peak of scores'
                        Ihat(:,side) = ea_discfibers_getpeak(v.*efieldT, obj.posvisible, obj.negvisible, 'peak');
                    case 'peak 5% of scores'
                        Ihat(:,side) = ea_discfibers_getpeak(v.*efieldT, obj.posvisible, obj.negvisible, 'peak5');
                end
            end
            Ihat = ea_nanmean(Ihat,2);
            [Rpred, Rpred_p] = corr(obj.responsevar(patsel), Ihat, 'type', corrType, 'rows', 'pairwise');
        end

        function capParpool(obj)
            % Ensures the current parallel pool (if any) has no more than
            % obj.permtestMaxWorkers workers before a parfor loop starts.
            % Each worker holds its own full copy of obj (including
            % obj.results.efield), so an uncapped pool (one worker per core)
            % can exhaust system memory on large analyses.
            p = gcp('nocreate');
            if isempty(p)
                parpool(obj.permtestMaxWorkers);
            elseif p.NumWorkers > obj.permtestMaxWorkers
                delete(p);
                parpool(obj.permtestMaxWorkers);
            end
        end
    end

    methods (Access = private,  Static)
        function vals = maskvals(vals, posvisible, negvisible)
            if ~posvisible
                vals(vals>0) = nan;
            end
            if ~negvisible
                vals(vals<0) = nan;
            end
        end
    end
end
