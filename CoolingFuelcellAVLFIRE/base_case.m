function [CH4, H2O, air] = base_case(T_env, H_env)
CH4 = struct(); CH4.phase = "gas"; CH4.n = SLM_to_mol_sec(3.954); CH4.yCH4 = 1; CH4.T = T_env;
H2O = struct(); H2O.phase = "liq"; H2O.yH2Ol = 1; H2O.T = T_env; H2O.n = CH4.n*2;
air = struct(); air.phase = "gas"; air.n = CH4.n*35; air.T = T_env;
air.yH2O = H_env*cal_yH2Osat(T_env); air.yO2 = 0.21*(1-air.yH2O); air.yN2 = 0.79*(1-air.yH2O);
end