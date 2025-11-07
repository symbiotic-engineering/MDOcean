function [analysis_list] = get_all_analyses()
%GET_ALL_ANALYSES Returns cell array of strings of all analysis classes
    dont_include = ["GenericAnalysis.m","analysis_to_calkit.m","get_all_analyses.m"];
    analysis_list = dir('analysis/*.m');
    if isempty(analysis_list)
        analysis_list = dir('mdocean/analysis/*.m');
        if isempty(analysis_list)
            analysis_list = dir('../mdocean/analysis/*.m');
            if isempty(analysis_list)
                msg = ['No analyses found. Current dir: ' pwd];
                ME = MException('mdocean:get_all_analyses:no_analyses', msg);
                throw(ME);
            end
        end
    end
    analysis_list = {analysis_list.name};
    analysis_list = analysis_list(~contains(analysis_list, dont_include));
end
