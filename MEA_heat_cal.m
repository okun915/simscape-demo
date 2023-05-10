function stack_T = MEA_heat_cal(heat_produced,pip_Q)

stack_mea_cp =  870; % [J/(kg*K)] Overall specific heat of membrane electrode assembly
stack_T = (heat_produced-pip_Q)/stack_mea_cp;

end
