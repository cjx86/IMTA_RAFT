function[f_p, S_p, conf, S_area, DOF] = spectra_cs_han(ts,del,n_ens,n_bands,alph)
%% Obtaining the Spectra of a time series
%% This Funtion is able to ensemble average and band average
%% ts is the time series input
%% del is the time spacing between samples where  f = 1/del

%% Ensemble Average

N = length(ts);
L_ens = floor(N/n_ens);

DOF = n_ens*n_bands*2;

% Break up the T.S. into the correct sized ensembles

j = 1;
ens = zeros(n_ens,L_ens);


check = hanning(L_ens);
% length(check)
% N
% pause

for i=1:L_ens:n_ens*L_ens
    
   ens(j,:)=ts(i:j*L_ens);
   ens_han(j,:) = ts(i:j*L_ens).*check';
   ens_var(j,:) = var(ens(j,:));
   ens_var_han(j,:) = var(ens_han(j,:));
   
   j = j+1;
   
end

% Remove the mean and take the F.T. of each of the ensembles

    % preallocate for speed
G_ens = zeros(n_ens,L_ens);
G_ensc = zeros(n_ens,L_ens);
ens_mn = zeros(n_ens,L_ens);
%f_j = zeros(n_ens,floor(L_ens/2-1));
S_j = zeros(n_ens,floor(L_ens/2-1));

for ii = 1:n_ens

[ens_mn(ii,1),~,~] = stats(ens(ii,:));    
ens(ii,:) = ens(ii,:) - ens_mn(ii,1);


    % F.T. of the data
G_ens(ii,:) = 1/L_ens*fft(ens(ii,:));

G_ensc(ii,:) = conj(G_ens(ii,:));

G_en(ii,:) = G_ens(ii,2:floor(L_ens/2));
G_enc(ii,:) = G_ensc(ii,2:floor(L_ens/2));

    % Spectrum of the data
S_j(ii,:) = 2*L_ens*del*(G_en(ii,:).*G_enc(ii,:)).*(ens_var./ens_var_han);


end

f_j = (1:L_ens/2)/(L_ens*del);

%Sum up the Spectrum of each ensemble
Sj_tot = zeros(1,floor(L_ens/2-1));
f_ens = f_j(1,:);

for p = 1:floor(L_ens/2-1)

    for jj = 1:n_ens
    
        Sj_tot(1,p) = Sj_tot(1,p) + S_j(jj,p);
       
    end

end

Sj_ens = Sj_tot/n_ens; 


%% Band Average The ensemble averaged Spectrum
clear i
clear j

L_b = floor(length(Sj_ens)/n_bands);
j = 1;

band = zeros(L_b,n_bands);

%Fill the array of bands up 
for i=1:n_bands:(n_bands*L_b)
   
   band(j,:)=Sj_ens(i:j*n_bands);
   j = j+1;
    
end

Sj_band = zeros(1,L_b);
clear k

%Sum up each band array
for k = 1:L_b
    
    for kk = 1:n_bands
        Sj_band(1,k) = Sj_band(1,k) + band(k,kk);
    end

end

%Create the frequency array

clear j
j = 1:L_b;
m1 = 1+(j-1)*n_bands;
m2 = j*n_bands;

f_p = (f_ens(m2)+f_ens(m1))./2;
    
S_p = Sj_band/n_bands; 

f_p = f_p(2:end);
S_p = S_p(2:end);

S_area=trapz(f_p,S_p);

[chi2l, chi2h]=chi2_vals(alph,DOF);
Sl=DOF*S_p/chi2l; %Random S^ value
Sh=DOF*S_p/chi2h;
conf=(Sh-Sl)';
    
end
