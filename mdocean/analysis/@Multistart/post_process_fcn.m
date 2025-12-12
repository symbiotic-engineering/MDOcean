function [fig_array,...
                 tab_array_display,...
                 tab_array_latex,...
                 end_result_struct,...
                 tab_firstrows,...
                 tab_colspecs] = post_process_fcn(intermed_result_struct)
            
            p = intermed_result_struct.p;
            b = intermed_result_struct.b;
            X_opt =intermed_result_struct.X_opt;
            objs = intermed_result_struct.objs;
            flags = intermed_result_struct.flags;
            x0s = intermed_result_struct.x0s;
            num_runs = intermed_result_struct.num_runs;

            [treeFig, parallelFig, barFig, results] = multistart_postpro(p,b,X_opt,objs,flags,x0s,num_runs);

            fig_array = [treeFig, parallelFig, barFig];
            
            tab_array_display = {results};
            tab_array_latex = {results};

            tab_firstrows = {[]};
            tab_colspecs = {[]};

            end_result_struct = struct();
        end