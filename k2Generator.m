function Q = k2Generator(N,f,spar)
if spar == 0
    Q(N+1,N+1) = 0;
    for i = 2:N
        Q(i,i+1) = (i-1)/N*(f/((N-i+1)/N + f*((i-1)/N)))*(N-i+1)/N;
        Q(i,i-1) = (N-i+1)/N*(1/((N-i+1)/N + f*((i-1)/N)))*(i-1)/N;
        Q(i,i) = -Q(i,i+1) - Q(i,i-1);
    end
    Qs = Q(2:N,2:N);
    V = Q(2:N,[1,N+1]);
    e1 = [1,zeros(1,N-2)];
    I = eye(N-1,N-1);
else
    Norig = N;
    N = floor(N/6);
    parfor k = 1:6
        vvdown = zeros(N-1,1);
        vvup = vvdown;
        vvstay = vvdown;
        count = 0;
        minrange = N*(k-1)+1;
        if k < 6
        maxrange = N*k;
        else
            maxrange = Norig;
        end
        %ii(1,3*(N-1)) = 0;
        %jj(1,3*(N-1)) = 0;
        %vv(1,3*(N-1)) = 0;
        for i = minrange:maxrange
            count = count + 1;
            vvdown(count) = (Norig-i)/Norig*(1/((Norig-i)/Norig + f*((i)/Norig)))*(i)/Norig;
            vvup(count) = (i)/Norig*(f/((Norig-i)/Norig + f*((i)/Norig)))*(Norig-i)/Norig;
            vvstay(count) = -vvdown(count)-vvup(count);
        end
        vvp{k} = [vvup,vvstay,vvdown];
        %Q{k} = sparse(ii,jj,vv);
    end
    vv = [0,0,0];
    for k = 1:6
        vv = [vv;vvp{k}];
    end
    Q = spdiags(vv,-1:1,Norig+1,Norig+1)';
end
