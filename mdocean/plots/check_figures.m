[fig_success,tab_success, fig_output, tab_output,fig_runtime,tab_runtime] = all_figures();

fig_success_copy = fig_success;
idx_success = cellfun(@isempty,fig_success);
fig_success_copy(idx_success) = {MException('','ok')};

msg = cellfun(@(c) c.message, fig_success_copy, 'UniformOutput', false);
stack = cellfun(@(c) c.stack, fig_success_copy, 'UniformOutput', false);

file = cell(size(fig_success_copy));
top_stack_cell = cellfun(@(x) x(1), stack(~idx_success), 'UniformOutput', false);
top_stack = [top_stack_cell{:}];
file(~idx_success) = {top_stack.name};

line = zeros(size(fig_success_copy));
line(~idx_success) = [top_stack.line];
table((1:length(msg)).', msg, file, line, 'VariableNames', {'#','Error Message','File','Line'})
