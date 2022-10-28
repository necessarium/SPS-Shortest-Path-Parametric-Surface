function SPS(S, u, v, ui, uf, vi, vf)
    %surface - points - plot

    [PI, PF] = deal(S(ui, vi),  S(uf, vf));
    mk_plt_pckg = @(x,y,z){x,y,z,'.','MarkerSize',20};
    PI_pck = mk_plt_pckg(PI(1), PI(2), PI(3));
    plot3(PI_pck{:}); hold on; stem3(PI_pck{1:3});
    PF_pck = mk_plt_pckg(PF(1), PF(2), PF(3));
    plot3(PF_pck{:}); hold on; stem3(PF_pck{1:3});

    a = ui-1; b = uf+1; a_b = abs(ui+uf)/10; 
    c = vi-1; d = vf+1; c_d = abs(vi+vf)/10;
    [U, V] = meshgrid(a:a_b:b, c:c_d:d);
    C_Sur = S(U, V); for i = 1:length(C_Sur)
    C_Sur{i} = double(C_Sur{i}); end
    [X, Y, Z] = C_Sur{:};
    surfc(X, Y, Z); hold on
    styl = {'LineWidth', 2};

    function SPS_Discrete
        %pre-loop setup
        r0 = PI; h0 = 0.23; x0 = [0, 0];
        rn = r0; hn = h0; xn = x0; [un, vn] = deal(ui, vi);
        mk_quiv_pckg = @(a, b)[a, b, {'off', 'r'}, styl];
        Su(u,v) = diff(S,u); Sv(u,v) = diff(S,v);
        
        while norm(rn-PF) > hn/2
            axis equal;
            an = sym('an', 'real');
            bn = sym('bn', 'real');
            kn = sym('kn', 'real');
            Sn_u = double(Su(un,vn)); 
            Sn_v = double(Sv(un,vn));
            Nn = cross(Sn_u, Sn_v);
        
            %find T
            Tn = an*Sn_u + bn*Sn_v; M_Tn = matlabFunction(Tn)
            sp = dot(Tn, PF-rn); M_sp = matlabFunction(sp); 
            O_sp = @(x)-M_sp(x(1), x(2));
            nonlcon = @(x)deal([], double(norm(M_Tn(x(1), x(2)))^2 - 1));
            xn = fmincon(O_sp, xn, [], [], [], [], [], [], nonlcon);
            [an, bn] = deal(xn(1), xn(2)); Tn = an*Sn_u + bn*Sn_v;
        
            %plot rn
            rn_plt_pckg = mk_plt_pckg(rn(1), rn(2), rn(3));
            plot3(rn_plt_pckg{:}); hold on
            
            %T step - plot
            C_rn = num2cell(rn); C_hTn = num2cell(hn*Tn);
            T_stp_pckg = mk_quiv_pckg(C_rn, C_hTn);
            quiver3(T_stp_pckg{:}); hold on
            
            %r at n+1 - solve next un, vn, kn
            rnp1 = rn + hn*Tn + kn*Nn;
            M_rnp1 = matlabFunction(rnp1);
            F = @(x)double(M_rnp1(x(1)) - S(x(2), x(3)));
            C_x_kuvn = num2cell(fsolve(F, [0, un, vn]));
            [kn, un, vn] = C_x_kuvn{:};
        
            %N step - plot
            rphTn = rn + hn*Tn;
            C_rphTn = num2cell(rphTn); C_kNn = num2cell(kn*Nn);
            N_stp_pckg = mk_quiv_pckg(C_rphTn, C_kNn);
            quiver3(N_stp_pckg{:});

            %update rn to rnp1
            rn = double(subs(rn + hn*Tn + kn*Nn, 'kn', kn));
        end
    end

    function SPS_Continuous(m, k)
        t = sym('t', 'real');
        T = (t.^((0:m).'));
        
        function I = min_int(x)  
            C = reshape(x, [2 m+1]);
            r = C*T; C_r = num2cell(r);
            Sr_p = diff(S(C_r{:}), t);
            m_p = simplify(norm(Sr_p)^k);
            M_intgd = matlabFunction(m_p);
            I = integral(M_intgd, 0, 1);
        end 

        function L = len_int(x, a, b)  
            C = reshape(x, [2 m+1]);
            r = C*T; C_r = num2cell(r);
            Sr_p = diff(S(C_r{:}), t);
            l_p = simplify(norm(Sr_p));
            M_intgd = matlabFunction(l_p);
            L = integral(M_intgd, a, b);
        end
        
        function [c,ceq] = varcon(x)
            C = reshape(x, [2 m+1]);
            r = matlabFunction(C*T);
            ceq1 = r(0).' - [ui vi];
            ceq2 = r(1).' - [uf vf];
            ceq = [ceq1 ceq2]; c = [];
        end
        
        C0 = zeros(2, m+1);
        C0(:,1:2) = [ui uf-ui; vi vf-vi];
        x0 = reshape(C0, [1 2*(m+1)]); 
        con_pckg = {[],[],[],[],[],[], @varcon};
        x_min = fmincon(@min_int,x0,con_pckg{:});

        len_int(x_min, 0, 1)
        C_min = reshape(x_min, [2 m+1]);

        R(t) = C_min*T; 
        Ce_rmin = num2cell(R(t));
        Ce_Sr = num2cell(S(Ce_rmin{:}));
        fplt_pckg = [Ce_Sr, {[0 1]}, styl];
        axis equal; fplot3(fplt_pckg{:})
    end

    SPS_Continuous(2, 2);
    SPS_Discrete();
end