function [diagnostic] = save_fig_with_diagnostic(fig, fig_name, pdf_prefix)
% saves a pdf of a figure, robustly handling cases of gobjects, deleted
% handles, and invalid figures

    if nargin<3
        pdf_prefix = '';
    end
    pdf_name = fullfile(pdf_prefix, fig_name);

    if isgraphics(fig) % if figure exists (didn't error first and wasn't deleted)
        if ~isempty(fig.UserData)
            % pdf already exists in files, just copy to folder
            copyfile(fig.UserData, pdf_name)
        else
            % save pdf from matlab figure output
            save_pdf(fig,pdf_name)
        end

        if nargout>0
            % in either case, use figure itself, not pdf, for printing the diagnostic
            diagnostic = matlab.unittest.diagnostics.FigureDiagnostic(fig,'Prefix',fig_name + '_');
        end

    elseif ~isvalid(fig)
        msg = 'Figure has been deleted, probably because of a "close all" in a subsequent script.';
        err = MException('MDOcean:deletedFigure',msg);
        throw(err)

    else % placeholder gobjects figure - don't save
        if nargout > 0
            diagnostic = matlab.unittest.diagnostics.DisplayDiagnostic('');
        end
    end
end
