prediction = @(V_hull) V_hull*.27876 + [-1,0,1]*280;

V_hull_nom = 1264; % V_f_tot from breakpoint in geometry
m_struct_nom = 213 % mass_f (metric tons)
m_struct_pred_nom  = prediction(V_hull_nom)

V_hull_opt = 7512;
m_struct_opt = 80
m_struct_pred_opt = prediction(V_hull_opt)
