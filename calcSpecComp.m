% predict for a multicomponent alloy system

clear;

load('EleIndex.mat');

% Fig 9
data1 = xlsread('Multi_BMG.xlsx','MGComp','C2:J352');  
data5 = xlsread('Multi_BMG.xlsx','MGEle','B2:I352'); 
% ------------------

% critical size
% data1 = xlsread('Multi_BMG.xlsx','GFA','H2:J231');  
% data5 = xlsread('Multi_BMG.xlsx','GFA','E2:G231');
% data6 = xlsread('Multi_BMG.xlsx','GFA','L2:L231');
% -------------------

% revised manuscript Fig 9
% data1 = xlsread('Multi_BMG.xlsx','CryComp','C1:H35');  
% data5 = xlsread('Multi_BMG.xlsx','CryEle','B1:G35');
% -------------------

nData = length(data1);
comp = cell(nData,1);
ele  = cell(length(data5),1);

% load properties
data2 = xlsread('properties20190429.xlsx','B2:E59');   
data3 = xlsread('mixingenthalpy20190328.xlsx','B2:BG59'); 
data4 = xlsread('properties20190429.xlsx','F2:N59'); 

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
zeta = 1/(1-xi);

% initialization
ME = zeros(nData,1);
delta = zeros(nData,1);
density = zeros(nData,1);
Sc_n = zeros(nData,1);
S_s_n = zeros(nData,1);
volume = zeros(nData,1);
mod_E = zeros(nData,1);
mod_G = zeros(nData,1);
mod_K = zeros(nData,1);
Tx = zeros(nData,1);
Tg = zeros(nData,1);
Tl = zeros(nData,1);
dTx = zeros(nData,1);
dchi = zeros(nData,1);
VEC = zeros(nData,1);
m = zeros(nData,1);
gamma_m = zeros(nData,1);
r_ave = zeros(nData,1);
tmp1 = zeros(nData,1);
w_total = zeros(nData,1);
w = zeros(nData,1);
sigma2 = zeros(nData,1);
sigma3 = zeros(nData,1);
tmp3 = zeros(nData,1);
tmp4 = zeros(nData,1);
y1 = zeros(nData,1);
y2 = zeros(nData,1);
y3 = zeros(nData,1);
tmp5 = zeros(nData,1);
mod_E_up = zeros(nData,1);
tmp6 = zeros(nData,1);
mod_E_lo = zeros(nData,1);
tmp7 = zeros(nData,1);
mod_K_up = zeros(nData,1);
tmp8 = zeros(nData,1);
mod_K_lo = zeros(nData,1);
tmp9 = zeros(nData,1);
mod_G_up = zeros(nData,1);
tmp10 = zeros(nData,1);
mod_G_lo = zeros(nData,1);
chi_ave = zeros(nData,1);
tmp11 = zeros(nData,1);

for kk=1:length(data1)
    
    comp{kk} = rmmissing(data1(kk,:));
    ele{kk} = rmmissing(data5(kk,:));
    
    if abs(sum(comp{kk})-100)>1e-5
%         error('composition sum should be 100');
        dchi(kk) = NaN;
        delta(kk) = NaN;
        density(kk) = NaN;
        dTx(kk) = NaN;
        gamma_m(kk) = NaN;
        m(kk) = NaN;
        ME(kk) = NaN;
        S_s_n(kk) = NaN;
        Sc_n(kk) = NaN;
        mod_K(kk) = NaN;
        VEC(kk) = NaN;
        continue
    else
        comp{kk} = comp{kk}/100;
    end
    if length(ele{kk})~=length(comp{kk})
%         error('number of composition and element should be equal');
        dchi(kk) = NaN;
        delta(kk) = NaN;
        density(kk) = NaN;
        dTx(kk) = NaN;
        gamma_m(kk) = NaN;
        m(kk) = NaN;
        ME(kk) = NaN;
        S_s_n(kk) = NaN;
        Sc_n(kk) = NaN;
        mod_K(kk) = NaN;
        VEC(kk) = NaN;
        continue
    end
    
    for ii=1:length(ele{kk})
        
        for jj=ii+1:length(ele{kk})
            ME(kk) = ME(kk) + 4*comp{kk}(ii)*comp{kk}(jj)*mixH(ele{kk}(ii),ele{kk}(jj));
        end

        r_ave(kk) = r_ave(kk) + comp{kk}(ii)*radius(ele{kk}(ii));
        if ii==length(ele{kk})
            for jj=1:length(ele{kk})
                tmp1(kk) = tmp1(kk)+ (1-radius(ele{kk}(jj))/r_ave(kk))^2*comp{kk}(jj);
            end
            delta(kk) = sqrt(tmp1(kk));
        end

        w_total(kk) = w_total(kk) + comp{kk}(ii)*a_Mass(ele{kk}(ii));
        if ii==length(ele{kk})
            for jj=1:length(ele{kk})
                w(kk) = w(kk) + comp{kk}(jj)*a_Mass(ele{kk}(jj))/w_total(kk)/atomicDensity(ele{kk}(jj)); % wt.% of i-th element
            end
            density(kk) = 1/w(kk)/1000;
        end

        Sc_n(kk) = Sc_n(kk) - fillmissing(comp{kk}(ii)*log(comp{kk}(ii)),'constant',0);

        sigma2(kk) = sigma2(kk) + comp{kk}(ii)*covDia(ele{kk}(ii))^2;
        sigma3(kk) = sigma3(kk) + comp{kk}(ii)*covDia(ele{kk}(ii))^3;
        for jj=ii+1:length(ele{kk})
            tmp3(kk) = tmp3(kk) + (covDia(ele{kk}(ii))+covDia(ele{kk}(jj)))*(covDia(ele{kk}(ii))-covDia(ele{kk}(jj)))^2*comp{kk}(ii)*comp{kk}(jj); 
            tmp4(kk) = tmp4(kk) + covDia(ele{kk}(ii))*covDia(ele{kk}(jj))*(covDia(ele{kk}(ii))-covDia(ele{kk}(jj)))^2*comp{kk}(ii)*comp{kk}(jj);
        end
        y1(kk) = 1/sigma3(kk)*tmp3(kk);
        y2(kk) = sigma2(kk)/sigma3(kk)^2*tmp4(kk);
        y3(kk) = sigma2(kk)^3/sigma3(kk)^2;
        S_s_n(kk) = (1.5*(zeta^2-1)*y1(kk)+1.5*(zeta-1)^2*y2(kk)-(0.5*(zeta-1)*(zeta-3)+log(zeta))*(1-y3(kk)));

        volume(kk) = volume(kk) + comp{kk}(ii)*a_Mass(ele{kk}(ii))/atomicDensity(ele{kk}(ii))*1000;
        tmp5(kk) = tmp5(kk) + comp{kk}(ii)*ele_mod_E(ele{kk}(ii))*a_Vol(ele{kk}(ii));
        mod_E_up(kk) = tmp5(kk)/volume(kk);
        tmp6(kk) = tmp6(kk) + comp{kk}(ii)*a_Vol(ele{kk}(ii))/ele_mod_E(ele{kk}(ii));
        mod_E_lo(kk) = volume(kk)/tmp6(kk);
        mod_E(kk) = (mod_E_up(kk) + mod_E_lo(kk))/2;

        tmp7(kk) = tmp7(kk) + comp{kk}(ii)*ele_mod_K(ele{kk}(ii))*a_Vol(ele{kk}(ii));
        mod_K_up(kk) = tmp7(kk)/volume(kk);
        tmp8(kk) = tmp8(kk) + comp{kk}(ii)*a_Vol(ele{kk}(ii))/ele_mod_K(ele{kk}(ii));
        mod_K_lo(kk) = volume(kk)/tmp8(kk);
        mod_K(kk) = (mod_K_up(kk) + mod_K_lo(kk))/2;

        tmp9(kk) = tmp9(kk) + comp{kk}(ii)*ele_mod_G(ele{kk}(ii))*a_Vol(ele{kk}(ii));
        mod_G_up(kk) = tmp9(kk)/volume(kk);
        tmp10(kk) = tmp10(kk) + comp{kk}(ii)*a_Vol(ele{kk}(ii))/ele_mod_G(ele{kk}(ii));
        mod_G_lo(kk) = volume(kk)/tmp10(kk);
        mod_G(kk) = (mod_G_up(kk) + mod_G_lo(kk))/2;

        Tx(kk) = 2.68*mod_K(kk) + 300;
        Tg(kk) = mod_E(kk)*2.5 + 220; 
        Tl(kk) = mod_E(kk)*volume(kk)/97.9/R*1000;
        dTx(kk) = Tx(kk) - Tg(kk);
        m(kk) = mod_K(kk)/mod_G(kk)*12 + 8;
        gamma_m(kk) = (2*Tx(kk)-Tg(kk))/Tl(kk);
        VEC(kk) = VEC(kk) + comp{kk}(ii)*a_VEC(ele{kk}(ii));

        chi_ave(kk) = chi_ave(kk) + comp{kk}(ii)*a_eNega(ele{kk}(ii));
        if ii==length(ele{kk})
            for jj=1:length(ele{kk})
                tmp11(kk) = tmp11(kk) + comp{kk}(jj)*(a_eNega(ele{kk}(jj))-chi_ave(kk))^2;
            end
            dchi(kk) = sqrt(tmp11(kk));
        end
    end
end
    
[x1,ic1] = normx(dchi',0,0.7,-1,1);
[x2,ic2] = normx(delta',0,0.2,-1,1);
[x3,ic3] = normx(density',3,15,-1,1);
[x4,ic4] = normx(dTx',-50,250,-1,1);
[x5,ic5] = normx(gamma_m',0.2,1.2,-1,1);
[x6,ic6] = normx(m',25,65,-1,1);
[x7,ic7] = normx(ME',-70,50,-1,1);
[x8,ic8] = normx(S_s_n',0,1,-1,1);
[x9,ic9] = normx(Sc_n',0,1.1,-1,1);
[x10,ic10] = normx(mod_K',10,200,-1,1); 
[x11,ic11] = normx(VEC',1,12,-1,1);

input = [x1' x2' x3' x4' x5' x6' x7' x8' x9' x10' x11']';% 

output = NN93(input)';

% critical size
% [GFA,tr] = rmmissing(data6);
% output = output(~tr);
% -----------------
