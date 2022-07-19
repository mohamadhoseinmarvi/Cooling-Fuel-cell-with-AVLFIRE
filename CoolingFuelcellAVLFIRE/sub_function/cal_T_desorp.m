function [T_desorp, ex_out] = cal_T_desorp(ex_in, T_desorp)
% Solve for the minimum desorption temperature and the condition of the exhaust leaving desorption stage.
% If exhaust leaving desorption yH2O is greater than the saturated water mole fraction, increase desorption temperature.
yH2O_ex_out = solve_yH2O_ex_out(ex_in, T_desorp);
while yH2O_ex_out > cal_yH2Osat(T_desorp)
    T_desorp = T_desorp+5;
    yH2O_ex_out = solve_yH2O_ex_out(ex_in, T_desorp);
end
[~, ex_out]=equation_set(yH2O_ex_out, ex_in, T_desorp);
end

function yH2O_ex_out = solve_yH2O_ex_out(ex_in, T_desorp)
% Unknowns: 1: ex_out.yH2O/y_EQ
options = optimset("MaxFunEvals", 10000, "MaxIter", 10000, 'Display', 'off');
% Initial guess:
y0 = ex_in.yH2O*1.1;
% Solve:
f = @(y) equation_set(y, ex_in, T_desorp);
[yH2O_ex_out, fval] = fsolve(f,y0,options);
if max(abs(fval))>1e-6
    return
end
end

function [F, ex_out]=equation_set(y, ex_in, T_desorp)
H2O_MW = 18.01528;  % [g/mol]
Q_adsorption = 2800;    

ex_out = cal_ex_out(ex_in, y, T_desorp);
dH = cal_stream_enthalpy(only_air(ex_in))-cal_stream_enthalpy(only_air(ex_out));
dmH2O_desorption = (ex_in.n*ex_in.yH2O-ex_out.n*ex_out.yH2O)*H2O_MW*Q_adsorption;
F = dH + dmH2O_desorption;
end

function ex_out = cal_ex_out(ex_in, yH2O_ex_out, T_desorp)
ex_out = ex_in;
ex_out.T = T_desorp;
ex_out.yH2O = yH2O_ex_out;
ex_out.n = ex_in.n*( 1-ex_in.yH2O )/( 1-ex_out.yH2O );
gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'};
for i = 1:length(gas_list)
    gas = gas_list{i};
    if isfield(ex_in, gas)
        ex_out.(gas) = ex_in.(gas)*ex_in.n/ex_out.n;
    end
end  
end

function out = only_air(in)
    out = in;
    out.n = in.n*(1-in.yH2O);
    gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(in, gas)
            out.(gas) = in.(gas)*in.n/out.n;
        end
    end
    out.yH2O = 0;
end