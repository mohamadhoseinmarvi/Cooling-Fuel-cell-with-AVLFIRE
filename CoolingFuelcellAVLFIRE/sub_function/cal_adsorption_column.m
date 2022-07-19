function [air_2, air_2_dynamic, T_bed_f, t_ads] = cal_adsorption_column(air_1, T_bed_0, Desiccant, yH2O_ads)
[~, liq_stream] = cal_isenthalpic_condensation(air_1);
if isstruct(liq_stream)
    fprintf("ERROR: input stream (in) is not all gaseous\n")
    return
end
if ~(isfield(air_1, "yH2O"))
    fprintf("ERROR: input stream (in) does not contain water")
end
H2O = struct();
H2O.Cp = 4.18; 
H2O.MW = 18.01528; 
mH2O_end = Desiccant.m*Desiccant.capacity;
tspan = [0 (Desiccant.m*Desiccant.capacity/(air_1.yH2O*air_1.n*H2O.MW*0.01))];
y0 = [0, T_bed_0];
opt = odeset('Events', @myEvent);
[t, y] = ode15s(@(t, y) find_prime(t, y, air_1), tspan, y0, opt);
t_ads = t(end);
T_bed_f = y(end,2);
T_bed_m = y(:,2);
air_2_dynamic = cal_out_stream(air_1, T_bed_m, yH2O_ads);
air_2_dynamic.t = t;
air_2 = cal_out_stream(air_1, T_bed_m, yH2O_ads);
air_2.T = trapz(t, y(:,2))/t_ads;

    function dy = find_prime(~, y, in)
    dy = zeros(size(y));
    T_bed = y(2);
    out = cal_out_stream(in, T_bed, yH2O_ads);
    dH = cal_stream_enthalpy(only_air(in))-cal_stream_enthalpy(only_air(out));
    dy(1) = ( in.n*in.yH2O - out.n*out.yH2O)*H2O.MW;
    dy(2) = ( dy(1)*Desiccant.Q_adsorption + dH )/( y(1)*H2O.Cp + Desiccant.m*Desiccant.Cp );
    end

    function [value, isterminal, direction] = myEvent(~, y)
    value      = y(1) - mH2O_end;
    isterminal = 1;  
    direction  = 0;
    end
end

function out = cal_out_stream(in, T_bed, yH2O_ads)
    out = struct();
    out.phase = "gas";
    out.yH2O = yH2O_ads;
    out.n = in.n*(1-in.yH2O)./(1-out.yH2O);
    gas_list = {'yCH4', 'yCO', 'yCO2', 'yH2', 'yO2', 'yN2'}; % yH2O not in list.
    for i = 1:length(gas_list)
        gas = gas_list{i};
        if isfield(in, gas)
            out.(gas) = in.(gas)*in.n./out.n;
        end
    end
    out.T = T_bed;
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