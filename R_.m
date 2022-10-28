function R_
    C_
     %[cos(u)*cos(v),sin(u)*cos(v), sin(v)]; [1, 5]; [-0.5, 0.5]
    %[u, v, v - sin(v)*u.^2 + v.^2]; [-0.75,-1.312]; [2.1,1.243];
    %[u, v, u.^2 + v.^2]; [-0.5, -1]; [1, 1.5];
    %[cos(u)+v,u+sin(v),v-sin(v)*u.^2+v.^2]; [-3.74,-2.32]; [2.11,4.7];
    %[v, u*sin(v), v*cos(u)]
    u = sym('u', 'real'); v = sym('v', 'real');
    S(u,v) = [cos(u)+v,u+sin(v),v-sin(v)*u.^2+v.^2];
    [ui, vi] = deal(-2.743, -2.32);
    [uf, vf] = deal(2.11, 1.7);
    SPS(S, u, v, ui, uf, vi, vf)
end