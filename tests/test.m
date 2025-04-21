classdef (SharedTestFixtures={ ...
        matlab.unittest.fixtures.CurrentFolderFixture('../mdocean')}) ...
        test < matlab.unittest.TestCase
    % class based unit tests, as in https://www.mathworks.com/help/matlab/matlab_prog/class-based-unit-tests.html
    
    properties (Constant)
        run_slow_tests = true;

        slow_figs = [16:25];
        slow_tabs = 1:8;
    end

    properties
        feasible_report
        failed_report
        simulated_report
        actual_report
        econ_fig_report

        feasible_wecsim
        failed_wecsim
        simulated_wecsim
        actual_wecsim
        econ_fig_wecsim

        uuid

        fig_success
        tab_success
        fig_output
        tab_output
        fig_runtime
        tab_runtime
    end

    % inputs for tests, including passing tolerances
    properties (TestParameter)
        field_report = fieldnames(validation_inputs('report'));
        field_wecsim = fieldnames(validation_inputs('wecsim'));
        rel_tol_report = {.1,.1,.1,.1,.01,.01,.25,.25,.25,.1,.1,.1,.1,.1,.1,.1,.1,.1};
        rel_tol_wecsim = {.01,.01,.01,.01, 0.1,0.1,.1,.1};
        which_figs = test.enumerateFigs()
        which_tabs = test.enumerateTabs()
    end

    % helper methods to enumerate all figures and tables
    methods (Static)
        function which_fig_struct = enumerateFigs()
            [~,~,~,~,~,~,num_figs,num_tabs,fig_names,~] = all_figures( [],[] );

            which_figs_vec = [1:num_figs zeros(1, num_tabs)];
            none = strcat(repmat({'none'},1,num_tabs), string(1:num_tabs));
            fig_names = matlab.lang.makeValidName([fig_names, none]);

            which_figs_cell = num2cell(which_figs_vec);
            which_fig_struct = cell2struct(which_figs_cell,fig_names,2);
        end
        function which_tab_struct = enumerateTabs()
            [~,~,~,~,~,~,num_figs,num_tabs,~,tab_names] = all_figures( [],[] );

            which_tabs_vec = [zeros(1,num_figs), 1:num_tabs];
            none = strcat(repmat({'none'},1,num_figs), string(1:num_figs));
            tab_names = matlab.lang.makeValidName([none, tab_names]);

            which_tabs_cell = num2cell(which_tabs_vec);
            which_tab_struct = cell2struct(which_tabs_cell,tab_names,2);
        end
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class

        function runNominalValidation(testCase)
            % this is a shared setup because the results are used by both
            % validateNominal and validateNominalFeasible
            [feas_r, fail_r, sim_r, act_r, ~, fig_r] = validate_nominal_RM3('report'); % report
            [feas_w, fail_w, sim_w, act_w, ~, fig_w] = validate_nominal_RM3('wecsim'); % wecsim
            
            testCase.feasible_report  = feas_r;
            testCase.failed_report    = fail_r;
            testCase.simulated_report = sim_r;
            testCase.actual_report    = act_r;
            testCase.econ_fig_report  = fig_r;

            testCase.feasible_wecsim  = feas_w;
            testCase.failed_wecsim    = fail_w;
            testCase.simulated_wecsim = sim_w;
            testCase.actual_wecsim    = act_w;
            testCase.econ_fig_wecsim  = fig_w;
        end

        function generateUUID(testCase)
            % generate unique identifier for each parallel worker
            % required to prevent file overalps for generated code
            testCase.uuid = parallel.pool.Constant(@() char(matlab.lang.internal.uuid()));
        end

        function runAllFigsTabs(testCase)
            [~,~,~,~,~,~,num_figs,num_tabs,...
                fig_names,tab_names] = all_figures( [], [] );

            all_figs = 1:num_figs;
            all_tabs = 1:num_tabs;
            if ~testCase.run_slow_tests
                run_figs = all_figs(~ismember(all_figs,testCase.slow_figs));
                run_tabs = all_tabs(~ismember(all_tabs,testCase.slow_tabs));
            else
                run_figs = all_figs;
                run_tabs = all_tabs;
            end

            [f_success_run,t_success_run,...
             f_output_run, t_output_run,...
             f_runtime_run,t_runtime_run] = all_figures( run_figs, run_tabs, testCase.uuid );

            % comment above and uncomment below to dry run
%             f_success_run = cell([1,length(run_figs)]);
%             t_success_run = cell([1,length(run_tabs)]);
%             f_output_run = gobjects(1,length(run_figs));
%             t_output_run = cell([1,length(run_tabs)]);
%             f_runtime_run = rand([1,length(run_figs)]);
%             t_runtime_run = rand([1,length(run_tabs)]);

            % store success info
            f_success = cell(1,num_figs);
            f_success(run_figs) = f_success_run;

            t_success = cell(1,num_tabs);
            t_success(run_tabs) = t_success_run;

            % store output
            f_output = gobjects(1, num_figs);
            f_output(run_figs) = f_output_run;

            t_output = cell(1,num_tabs);
            t_output(run_tabs) = t_output_run;

            % store runtimes
            f_runtime = NaN(1,num_figs);
            f_runtime(run_figs) = f_runtime_run;

            t_runtime = NaN(1,num_tabs);
            t_runtime(run_tabs) = t_runtime_run;

            % add runtime figure
            try
                times = [f_runtime t_runtime];
                names = remove_underscores([fig_names tab_names]);
                
                runtimeFig = figure;
                hold on
                idx_plot = false([1,length(times)]);
                idx_plot([15 16 22 26 end-1 end]) = true;
                names = reordercats(categorical(names(idx_plot)),names(idx_plot));
                bar(names,times(idx_plot));
                for i = 1:length(times)
                    if times(i) == 0 || ~isfinite(times(i))
                        plot(i,0,'rx')
                    end
                end
                title('Runtime (seconds)')
                hold off
                improvePlot
                set(runtimeFig,"Position",[1 41 1536 845])
                f_output(31) = runtimeFig;
            catch err
                f_success{31} = err;
            end


            % assign all in property
            testCase.fig_success = f_success;
            testCase.tab_success = t_success;
            testCase.fig_output  = f_output;
            testCase.tab_output  = t_output;
            testCase.fig_runtime = f_runtime;
            testCase.tab_runtime = t_runtime;

        end
    end

    methods(TestClassTeardown)
        function deleteGeneratedFiles(testCase)
            % Create a wildcard pattern
            pattern = fullfile('**',['*' testCase.uuid.Value '*']);
            
            % Get a list of folders that match the pattern
            matchingFiles = dir(pattern);
            foldersToDelete = matchingFiles([matchingFiles.isdir]);
            foldersToDelete = fullfile({foldersToDelete.folder},{foldersToDelete.name});

            % Delete matching folders
            for d = 1:length(foldersToDelete)
                rmpath(foldersToDelete{d})
                rmdir(foldersToDelete{d},'s');
            end

        end
    end
    
    % Test methods
    methods(Test, ParameterCombination='sequential')

        % run every figure and log it
        function allFiguresRun(testCase, which_figs, which_tabs)

            if ~test.run_slow_tests % mark slow tests as filtered
                testCase.assumeFalse( ismember(which_figs, test.slow_figs) );
                testCase.assumeFalse( ismember(which_tabs, test.slow_tabs) );
            end

            if which_figs ~= 0 % figure
                success_criterion = testCase.fig_success{which_figs};
                fig_out = testCase.fig_output(which_figs);

                fig_name = ['Figure_' num2str(which_figs)];
                pdf_name = ['../test-results/' fig_name];
                
                if isgraphics(fig_out) % if figure exists (didn't error first and wasn't deleted)
                    if ~isempty(fig_out.UserData)
                        % pdf already exists in files, just copy to folder
                        copyfile(fig_out.UserData, pdf_name)
                    else
                        % save pdf from matlab figure output
                        save_pdf(fig_out,pdf_name)
                    end
                    % in either case, use figure itself, not pdf, for printing the diagnostic
                    diagnostic = matlab.unittest.diagnostics.FigureDiagnostic(fig_out,'Prefix',[fig_name '_']);
                elseif ~isvalid(fig_out)
                    error('Figure has been deleted, probably because of a "close all" in a subsequent script.')
                end   

            else % table
                success_criterion = testCase.tab_success{which_tabs};
                tab_out = testCase.tab_output(which_tabs);

                diagnostic = matlab.unittest.diagnostics.DisplayDiagnostic(tab_out{:});
            end

            if iscell(success_criterion)
                success_criterion = success_criterion{:};
            end

            if isempty(success_criterion)
                success_criterion = 1;
            end

            if isa(success_criterion,'MException')
                if isempty(success_criterion.stack)
                    throw(success_criterion)
                else
                    rethrow(success_criterion)
                end
            else
                testCase.verifyGreaterThan(success_criterion, 0, diagnostic);
            end

        end

        function validateNominalReport(testCase, field_report, rel_tol_report)
            sim = testCase.simulated_report.(field_report);
            act = testCase.actual_report.(field_report);
            if strcmp(field_report,'LCOE')
                diagnostic = matlab.unittest.diagnostics.FigureDiagnostic(testCase.econ_fig_report,'Prefix','econ_validation_report');
            else
                diagnostic = '';
            end
            testCase.verifyEqual(sim, act, 'RelTol',rel_tol_report,diagnostic)
        end

        function validateNominalWecsim(testCase, field_wecsim, rel_tol_wecsim)
            sim = testCase.simulated_wecsim.(field_wecsim);
            act = testCase.actual_wecsim.(field_wecsim);
            if strcmp(field_wecsim,'LCOE')
                diagnostic = matlab.unittest.diagnostics.FigureDiagnostic(testCase.econ_fig_wecsim,'Prefix','econ_validation_wecsim');
            else
                diagnostic = '';
            end
            testCase.verifyEqual(sim, act, 'RelTol',rel_tol_wecsim,diagnostic)
        end

    end

    methods(Test)
        % Static tests
        function validateNominalReportFeasible(testCase)
            testCase.onFailure( ['Nominal design violates these constraints: ', testCase.failed_report] );
            testCase.verifyTrue(testCase.feasible_report);
        end

        function validateNominalWecsimFeasible(testCase)
            testCase.onFailure( ['Nominal design violates these constraints: ', testCase.failed_wecsim] );
            testCase.verifyTrue(testCase.feasible_wecsim);
        end

        function validateNominalHydroCoeffs(testCase)
            mean_err = hydro_coeff_err(false);
            testCase.verifyLessThanOrEqual(mean_err, 0.10)
        end

        function hydrodynamicLimitObeyed(testCase)
            ratio = check_max_CW(testCase.uuid.Value);
            ratio(isnan(ratio)) = 0;
            testCase.verifyLessThanOrEqual( ratio, 1 );
        end

    end
    
end
