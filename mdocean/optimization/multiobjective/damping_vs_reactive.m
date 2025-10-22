function damping_vs_reactive(p,b)

    p.control_type = 'damping';
    pareto_search(p,b, b.filename_uuid)

    p.control_type = 'reactive';
    pareto_search(p,b, b.filename_uuid)

    pareto_curve_heuristics(true)

end