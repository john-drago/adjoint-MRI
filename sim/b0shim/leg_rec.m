function out=leg_rec(n, m, u)
% compute legendre polynomial values for dat
% using RRI's recursive relations

if m>0
    p_mm=(-1)^m*ffactorial(2*m-1)*(1-u.^2).^(m/2);
else
    p_mm=1;
end

if (n==m)
    out=p_mm;
else
    p_mm1=(2*m+1)*u.*p_mm;
    
    if (n==m+1)
        out=p_mm1;
    else
        % recursive calculation needed
        a=m+2;
        p_ma_2=p_mm;
        p_ma_1=p_mm1;
        
        while true
            p_ma=((2*a-1)*u.*p_ma_1-(a+m-1)*p_ma_2)/(a-m);
            if a==n
                break;
            end
            % prepare next iteration
            p_ma_2=p_ma_1;
            p_ma_1=p_ma;
            a=a+1;
        end
        
        out=p_ma;
    end
end
