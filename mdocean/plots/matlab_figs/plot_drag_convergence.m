function fig = plot_drag_convergence(X_f_guesses, X_s_guesses, phase_X_f_guesses, phase_X_s_guesses, ...
                                    iters, multibody, X_tol, phase_X_tol)
    x_iter = 1:iters;
    fig = figure;
    plot(x_iter,X_f_guesses,'DisplayName','|X_f|')
    hold on
    plot(x_iter,phase_X_f_guesses,'DisplayName','\angleX_f')
    if multibody
        plot(x_iter,X_s_guesses,'DisplayName','|X_s|')
        plot(x_iter,phase_X_s_guesses,'DisplayName','\angleX_s')
    end
    legend

    x_iter_tol = x_iter(2:end-1)+1;
    plot(x_iter_tol,X_f_guesses(2:end-1)+X_tol,'k--',x_iter_tol,phase_X_f_guesses(2:end-1)+phase_X_tol,'k--','HandleVisibility','off')
    plot(x_iter_tol,X_f_guesses(2:end-1)-X_tol,'k--',x_iter_tol,phase_X_f_guesses(2:end-1)-phase_X_tol,'k--','HandleVisibility','off')
    if multibody
        plot(x_iter_tol,X_s_guesses(2:end-1)+X_tol,'k--',x_iter_tol,phase_X_s_guesses(2:end-1)+phase_X_tol,'k--','HandleVisibility','off')
        plot(x_iter_tol,X_s_guesses(2:end-1)-X_tol,'k--',x_iter_tol,phase_X_s_guesses(2:end-1)-phase_X_tol,'k--','HandleVisibility','off')
    end
    xlabel('Iterations')
        title('Drag Nonlinearity Convergence')
    improvePlot
end