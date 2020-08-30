clear all;
clc;

load('full.mat');
data2 = xlsread('properties20190429.xlsx','B2:E59');   
data3 = xlsread('mixingenthalpy20190328.xlsx','B2:BG59');  % read the mixing enthalpy in the database
data4 = xlsread('properties20190429.xlsx','F2:N59');   % read density (Column I), atomic mass (Column J)and covalent radii (Column K)

n = length(full_c_A);
nSys = length(full_ele_A);

full_ME = [];
full_delta = [];
full_density = [];
full_Sc_n = [];
full_S_s_n = [];
full_volume = [];
full_mod_E = [];
full_mod_G = [];
full_mod_K = [];
full_Tx = [];
full_Tg = [];
full_Tl = [];
full_dTx = [];
full_Trg = [];
full_gamma = [];
full_omega = [];
full_dchi = [];
full_VEC = [];
full_m = [];
full_alpha = [];
full_beta = [];
full_gamma_m = [];
full_paraxi = [];

c_A = full_c_A./100;
c_B = full_c_B./100;
c_C = full_c_C./100;

mixH = data3; % mixing enthalpy ()
radius = data2(:,1); % atomic radius (A)
a_Mass = data4(:,5); % atomic mass (g)
atomicDensity = data4(:,4); % atomic density (kg/m^3)
covDia = data4(:,6)*2; % covalent diameter (Angstrom)
ele_mod_E = data4(:,1); % Young's modulus (GPa)
ele_mod_G = data4(:,8); % shear modulus (GPa)
ele_mod_K = data4(:,2); % bulk modulus (GPa)
a_Vol = data4(:,9); % molar volume (cm^3/mol)
R = 8.314; % gas constant,J/K/mol
a_eNega = data2(:,3); % Pauling electronegativity
a_VEC = data2(:,4);% VEC
xi = 0.64; % packing fraction

for ii=1:nSys
    
    ele_A = repmat(full_ele_A(ii,1),n,1); % repeat each element
    ele_B = repmat(full_ele_B(ii,1),n,1);
    ele_C = repmat(full_ele_C(ii,1),n,1);

    % calculate mixing anthalpy
    ME = 4.*c_A.*c_B.*diag(mixH(ele_A,ele_B).*eye(n))...
        +4.*c_A.*c_C.*diag(mixH(ele_A,ele_C).*eye(n))...
        +4.*c_C.*c_B.*diag(mixH(ele_C,ele_B).*eye(n));
    full_ME = [full_ME;ME];

    % calculate atomic size difference
    r_ave = c_A.*radius(ele_A) + c_B.*radius(ele_B) + c_C.*radius(ele_C);
    delta = sqrt((1-radius(ele_A)./r_ave).^2.*c_A + (1-radius(ele_B)./r_ave).^2.*c_B + (1-radius(ele_C)./r_ave).^2.*c_C);
    full_delta = [full_delta;delta];

    % calculate density
    w_total = c_A.*a_Mass(ele_A) + c_B.*a_Mass(ele_B) + c_C.*a_Mass(ele_C);
    w_A = c_A.*a_Mass(ele_A)./w_total; % wt.% of A element
    w_B = c_B.*a_Mass(ele_B)./w_total; % wt.% of B element
    w_C = c_C.*a_Mass(ele_C)./w_total; % wt.% of C element
    density = 1./(w_A./atomicDensity(ele_A) + w_B./atomicDensity(ele_B) + w_C./atomicDensity(ele_C))/1000; % g/cm^3
    full_density = [full_density;density];

    % calculate configurational entropy
    Sc_A = c_A.*log(c_A);
    Sc_A(c_A==0) = 0;
    Sc_B = c_B.*log(c_B);
    Sc_B(c_B==0) = 0;
    Sc_C = c_C.*log(c_C);
    Sc_C(c_C==0) = 0;
    Sc_n = -(Sc_A + Sc_B + Sc_C); % configurational entropy normalized by gas constant R
    full_Sc_n = [full_Sc_n;Sc_n];

    sigma2 = c_A.*covDia(ele_A).^2 + c_B.*covDia(ele_B).^2 + c_C.*covDia(ele_C).^2; % unit:pm^2
    sigma3 = c_A.*covDia(ele_A).^3 + c_B.*covDia(ele_B).^3 + c_C.*covDia(ele_C).^3; % unit:pm^3
    y1 = 1./sigma3.*((covDia(ele_A)+covDia(ele_B)).*(covDia(ele_A)-covDia(ele_B)).^2.*c_A.*c_B...
        +(covDia(ele_B)+covDia(ele_C)).*(covDia(ele_B)-covDia(ele_C)).^2.*c_B.*c_C...
        +(covDia(ele_A)+covDia(ele_C)).*(covDia(ele_A)-covDia(ele_C)).^2.*c_A.*c_C); % dimensionless
    y2 = sigma2./sigma3.^2.*((covDia(ele_A).*covDia(ele_B)).*(covDia(ele_A)-covDia(ele_B)).^2.*c_A.*c_B...
        +(covDia(ele_B).*covDia(ele_C)).*(covDia(ele_B)-covDia(ele_C)).^2.*c_B.*c_C...
        +(covDia(ele_A).*covDia(ele_C)).*(covDia(ele_A)-covDia(ele_C)).^2.*c_A.*c_C); % dimensionless
    y3 = sigma2.^3./sigma3.^2; % dimensionless
    zeta = 1/(1-xi);
    S_s_n = (1.5*(zeta^2-1).*y1+1.5*(zeta-1)^2.*y2-(0.5*(zeta-1)*(zeta-3)+log(zeta)).*(1-y3)); % S_sigma normalized by kB
    full_S_s_n = [full_S_s_n;S_s_n];

    volume = (c_A.*a_Mass(ele_A)./atomicDensity(ele_A)...
        + c_B.*a_Mass(ele_B)./atomicDensity(ele_B)...
        + c_C.*a_Mass(ele_C)./atomicDensity(ele_C))*1000; % cm^3
    full_volume = [full_volume;volume];

    mod_E_up = (c_A.*ele_mod_E(ele_A).*a_Vol(ele_A)...
           + c_B.*ele_mod_E(ele_B).*a_Vol(ele_B)...
           + c_C.*ele_mod_E(ele_C).*a_Vol(ele_C))./volume;

    mod_E_lo = volume./(c_A.*a_Vol(ele_A)./ele_mod_E(ele_A)...
              +c_B.*a_Vol(ele_B)./ele_mod_E(ele_B)...
              +c_C.*a_Vol(ele_C)./ele_mod_E(ele_C));
    mod_E = (mod_E_up + mod_E_lo)/2;
    full_mod_E = [full_mod_E;mod_E];

    mod_K_up = (c_A.*ele_mod_K(ele_A).*a_Vol(ele_A)...
              + c_B.*ele_mod_K(ele_B).*a_Vol(ele_B)...
              + c_C.*ele_mod_K(ele_C).*a_Vol(ele_C))./volume;

    mod_K_lo = volume./(c_A.*a_Vol(ele_A)./ele_mod_K(ele_A)...
              +c_B.*a_Vol(ele_B)./ele_mod_K(ele_B)...
              +c_C.*a_Vol(ele_C)./ele_mod_K(ele_C));
    mod_K = (mod_K_up + mod_K_lo)/2;
    full_mod_K = [full_mod_K;mod_K];

    Tx = 2.68 * mod_K + 300;
    full_Tx = [full_Tx;Tx];

    Tg = mod_E * 2.5 + 220;
    full_Tg = [full_Tg;Tg];

    Tl= mod_E.*volume./97.9./R.*1000;
    full_Tl = [full_Tl;Tl];

    dTx = Tx - Tg;
    full_dTx = [full_dTx;dTx];

    Trg = Tg./Tl;
    full_Trg = [full_Trg;Trg];

    gamma = Tx./(Tg+Tl);
    full_gamma = [full_gamma;gamma];

    omega = Tg./Tx - 2*Tg./(Tg+Tl);
    full_omega = [full_omega;omega];

    chi_ave = c_A.*a_eNega(ele_A) + c_B.*a_eNega(ele_B) + c_C.*a_eNega(ele_C);
    dchi = sqrt(c_A.*(a_eNega(ele_A)-chi_ave).^2+c_B.*(a_eNega(ele_B)-chi_ave).^2+c_C.*(a_eNega(ele_C)-chi_ave).^2);
    full_dchi = [full_dchi;dchi];

    VEC = c_A.*a_VEC(ele_A) + c_B.*a_VEC(ele_B) + c_C.*a_VEC(ele_C);
    full_VEC = [full_VEC;VEC];

    mod_G_up = (c_A.*ele_mod_G(ele_A).*a_Vol(ele_A)...
              + c_B.*ele_mod_G(ele_B).*a_Vol(ele_B)...
              + c_C.*ele_mod_G(ele_C).*a_Vol(ele_C))./volume;

    mod_G_lo = volume./(c_A.*a_Vol(ele_A)./ele_mod_G(ele_A)...
              +c_B.*a_Vol(ele_B)./ele_mod_G(ele_B)...
              +c_C.*a_Vol(ele_C)./ele_mod_G(ele_C));
    mod_G = (mod_G_up + mod_G_lo)/2;
    full_mod_G = [full_mod_G;mod_G];

    m = mod_K./mod_G.*12 + 8;
    full_m = [full_m;m];
    
    alpha = Tx./Tl; % 
    full_alpha = [full_alpha;alpha];
    
    beta = Tx./Tg + Tg./Tl; % 
    full_beta = [full_beta;beta];
    
    gamma_m = (2*Tx-Tg)./Tl; % 
    full_gamma_m = [full_gamma_m;gamma_m];
    
    paraxi = Tg./Tl + dTx./Tx; % 
    full_paraxi = [full_paraxi;paraxi];
end

c_A = repmat(c_A,nSys,1);
c_B = repmat(c_B,nSys,1);
c_C = repmat(c_C,nSys,1);
ele_A = kron(full_ele_A,ones(n,1)); % repeat each element
ele_B = kron(full_ele_B,ones(n,1));
ele_C = kron(full_ele_C,ones(n,1));

save fullFeatures.mat c_A c_B c_C full_dchi full_delta full_density full_dTx ele_A ele_B ele_C...
    full_gamma full_m full_ME full_mod_E full_mod_G full_mod_K full_S_s_n full_Sc_n full_Tg...
    full_Tl full_Trg full_Tx full_VEC full_volume full_omega full_alpha full_beta full_gamma_m full_paraxi;