function [air_1, air_2, ex_1, T_desorp, t_cycle] = cal_adsorp_desorp_cycle(ex_FC, air_1, Desiccant, T_desorp, yH2O_ads)
y0 = initial_guess(ex_FC, air_1, Desiccant, T_desorp, yH2O_ads);
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
f = @(y) equation_set(y, air_1, ex_FC, Desiccant, T_desorp, yH2O_ads);
[n_air_1, fval] = fsolve(f, y0, options);

if max(abs(fval))>1e-6
    return
end

[~, air_1, air_2, ex_1, T_desorp, t_cycle] = equation_set(n_air_1, air_1, ex_FC, Desiccant, T_desorp, yH2O_ads);

end

function [F, air_1, air_2, ex_1, T_desorp, t_ads] = equation_set(y, air_1, ex_FC, Desiccant, T_desorp, yH2O_ads)
air_1.n = y;
T_desorp = cal_T_desorp(ex_FC, T_desorp);
[air_2, ~, T_bed_f, t_ads] = cal_adsorption_column(air_1, T_desorp, Desiccant, yH2O_ads);
[ex_1, T_desorp, t_desorp] = cal_desorption_column(ex_FC, T_bed_f, Desiccant, T_desorp);
F = t_ads - t_desorp;
end

function y0 = initial_guess(ex_FC, air_1, Desiccant, T_desorp, yH2O_ads)
[~, ex_out] = cal_T_desorp(ex_FC, T_desorp);
H2O_MW = 18.01528;  
air_1.n = 1;
d_ads_H2O_per_mole = air_1.n*air_1.yH2O - air_1.n*(1-air_1.yH2O)./(1-yH2O_ads)*yH2O_ads;
dH = cal_stream_enthalpy(ex_FC)-cal_stream_enthalpy(ex_out);
d_des_H2O = dH/Desiccant.Q_adsorption/H2O_MW;
y0 = d_des_H2O/d_ads_H2O_per_mole;
end