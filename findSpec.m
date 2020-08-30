% make predictions for a specific ternary system

clear;

load('fullFeatures.mat');
load('EleIndex.mat');

[x1,ic1] = normx(full_dchi',0,0.7,-1,1);
[x2,ic2] = normx(full_delta',0,0.2,-1,1);
[x3,ic3] = normx(full_density',3,15,-1,1);
[x4,ic4] = normx(full_dTx',-50,250,-1,1);
[x5,ic5] = normx(full_gamma_m',0.2,1.2,-1,1);
[x6,ic6] = normx(full_m',25,65,-1,1);
[x7,ic7] = normx(full_ME',-70,50,-1,1);
[x8,ic8] = normx(full_S_s_n',0,1,-1,1);
[x9,ic9] = normx(full_Sc_n',0,1.1,-1,1);
[x10,ic10] = normx(full_mod_K',10,200,-1,1);
[x11,ic11] = normx(full_VEC',1,12,-1,1);

input = [x1' x2' x3' x4' x5' x6' x7' x8' x9' x10' x11']'; 

% id = ele_A==Mg & ele_B==Al & ele_C==Ti;
% id = ele_A==Ti & ele_B==Cu & ele_C==Zr;
% id = ele_A==Cu & ele_B==Nb & ele_C==Hf; nAA='Cu';nBB='Nb';nCC='Hf'; %Q2
% id = ele_A==Al & ele_B==Nb & ele_C==La; nAA='Al';nBB='Nb';nCC='La'; %Q2
id = ele_A==Ni & ele_B==Zr & ele_C==Nb; nAA='Ni';nBB='Zr';nCC='Nb'; %Q2
% id = ele_A==Al & ele_B==Cu & ele_C==La; nAA='Al';nBB='Cu';nCC='La'; %Q2

id = ele_A==Ti & ele_B==Fe & ele_C==Cu; nAA='Ti';nBB='Fe';nCC='Cu';
% id = ele_A==Ti & ele_B==Ni & ele_C==Zr; nAA='Ti';nBB='Ni';nCC='Zr';
% id = ele_A==Ti & ele_B==Co & ele_C==Zr; nAA='Ti';nBB='Co';nCC='Zr';% Sci. Ad.
% id = ele_A==Fe & ele_B==Co & ele_C==Zr; nAA='Fe';nBB='Co';nCC='Zr';% Sci. Ad.

inputs = [input(:,id)];
c_AA = c_A(id);
c_BB = c_B(id);
c_CC = c_C(id);
ele_AA = ele_A(id);
ele_BB = ele_B(id);
ele_CC = ele_C(id);

YY = NN93(inputs)';

% figure
% colormap(jet)
% [hg,htick,hcb]=tersurf(c_AA,c_BB,c_CC,YY);
% hlabels=terlabel(nAA,nBB,nCC);
% caxis([0 1]);