function[corr] =  cross_corr(minlag,maxlag,x,y)
%% minlag and max lag are the desired LAG numbers
%% A POSITIVE lag max correlation value means that 'X' leads 'Y'
%% A NEGATIVE lag max correlation value means that 'Y' leads 'X'



N = length(x);
corr = zeros(1,(maxlag-minlag)+1);
q = 0;

for k = minlag:maxlag

Sx = 0;
Sylag = 0;
Sxx = 0;
Syylag = 0;
Sxy = 0;
    
    Np = N-k;
    
    m = 0;
        
    for ii = 1:Np
        
        if isnan(x(ii)) ~= 1 && isnan(y(ii+k)) ~= 1
            
            Sx = Sx + x(ii);
            Sxx = Sxx + x(ii)*x(ii);
            Sylag = Sylag + y(ii+k);
            Syylag = Syylag + y(ii+k)*y(ii+k);
            Sxy =  Sxy + x(ii)*y(ii+k);
            m = m+1;
            
        end

        S_xy = Sxy/m;
        S_x = Sx/m;
        S_ylag = Sylag/m;
        S_xx = Sxx/m;
        S_yylag = Syylag/m;
    
    end
    
    q = q+1;
    
    corr(q) = (S_xy-S_x*S_ylag)/(sqrt(S_xx-S_x^2)*sqrt(S_yylag-S_ylag^2));
    
end

end