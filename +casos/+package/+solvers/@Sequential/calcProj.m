
function sol_proj = calcProj(obj,sos_g,idx)
    sol_proj = obj.projCon('p',sos_g(idx));
end