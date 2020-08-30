% Neural network model training

clear;

% load dataset
load('Dataset.mat');

% normalized input data
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

input = [x1' x2' x3' x4' x5' x6' x7' x8' x9' x10' x11']';
target = phase';

NNname = ['NN'];

% Create a BP Neural Network
hiddenLayerSize = 22; 

net = feedforwardnet(hiddenLayerSize,'trainbfg');

net.layers{1}.transferFcn = 'tansig'; % activation function
net.trainParam.max_fail = 10;
net.trainParam.goal = 0.05;
net.biases{1,1}.initFcn = 'rands';
net.biases{2,1}.initFcn = 'rands';
net.inputWeights{1,1}.initFcn = 'rands';
net.layerWeights{2,1}.initFcn = 'rands';

% Train the Network
[net,tr] = train(net,input,target);
genFunction(net,NNname);
output = net(input);

% Evaluate the model
[c,cm,ind,per]= confusion(target,output);
TN = cm(1,1);    FN = cm(2,1);
FP = cm(1,2);    TP = cm(2,2);

ACC = (TN+TP)/(TN+FN+FP+TP);

Pr_P = TP/(TP+FP);
Re_P = TP/(TP+FN);
F1_P = 2*Pr_P*Re_P/(Pr_P+Re_P);
Pr_N = TN/(TN+FN);
Re_N = TN/(TN+FP);
F1_N = 2*Pr_N*Re_N/(Pr_N+Re_N);

p_o = (TN+TP)/(TN+FN+FP+TP);
p_e = ((TN+FP)*(TN+FN)+(FN+TP)*(FP+TP))/(TN+FN+FP+TP)^2;
C_kappa = (p_o-p_e)/(1-p_e);

[prec, tpr, fpr, npv, thresh] = prec_rec(output, target);
% --------------

% test dataset
testTarget = target(tr.testInd);
testOutput = output(tr.testInd);
[cT,cmT,indT,perT]= confusion(testTarget,testOutput);

TNt = cmT(1,1);    FNt = cmT(2,1);
FPt = cmT(1,2);    TPt = cmT(2,2);

ACCt = (TNt+TPt)/(TNt+FNt+FPt+TPt);

Pr_Pt = TPt/(TPt+FPt);
Re_Pt = TPt/(TPt+FNt);
F1_Pt = 2*Pr_Pt*Re_Pt/(Pr_Pt+Re_Pt);
Pr_Nt = TNt/(TNt+FNt);
Re_Nt = TNt/(TNt+FPt);
F1_Nt = 2*Pr_Nt*Re_Nt/(Pr_Nt+Re_Nt);

p_ot = (TNt+TPt)/(TNt+FNt+FPt+TPt);
p_et = ((TNt+FPt)*(TNt+FNt)+(FNt+TPt)*(FPt+TPt))/(TNt+FNt+FPt+TPt)^2;
C_kappat = (p_ot-p_et)/(1-p_et);

[prect, tprt, fprt, npvt, thresht] = prec_rec(testOutput, testTarget);
% -----------

% validation dataset
valTarget = target(tr.valInd);
valOutput = output(tr.valInd);
[cV,cmV,indV,perV]= confusion(valTarget,valOutput);
TNv = cmV(1,1);    FNv = cmV(2,1);
FPv = cmV(1,2);    TPv = cmV(2,2);

ACCv = (TNv+TPv)/(TNv+FNv+FPv+TPv);

Pr_Pv = TPv/(TPv+FPv);
Re_Pv = TPv/(TPv+FNv);
F1_Pv = 2*Pr_Pv*Re_Pv/(Pr_Pv+Re_Pv);
Pr_Nv = TNv/(TNv+FNv);
Re_Nv = TNv/(TNv+FPv);
F1_Nv = 2*Pr_Nv*Re_Nv/(Pr_Nv+Re_Nv);

p_ov = (TNv+TPv)/(TNv+FNv+FPv+TPv);
p_ev = ((TNv+FPv)*(TNv+FNv)+(FNv+TPv)*(FPv+TPv))/(TNv+FNv+FPv+TPv)^2;
C_kappav = (p_ov-p_ev)/(1-p_ev);

[precv, tprv, fprv, npvv, threshv] = prec_rec(valOutput, valTarget);
% -----------

% train dataset
trainTarget = target(tr.trainInd);
trainOutput = output(tr.trainInd);
[cR,cmR,indR,perR]= confusion(trainTarget,trainOutput);
TNr = cmR(1,1);    FNr = cmR(2,1);
FPr = cmR(1,2);    TPr = cmR(2,2);

ACCr = (TNr+TPr)/(TNr+FNr+FPr+TPr);

Pr_Pr = TPr/(TPr+FPr);
Re_Pr = TPr/(TPr+FNr);
F1_Pr = 2*Pr_Pr*Re_Pr/(Pr_Pr+Re_Pr);
Pr_Nr = TNr/(TNr+FNr);
Re_Nr = TNr/(TNr+FPr);
F1_Nr = 2*Pr_Nr*Re_Nr/(Pr_Nr+Re_Nr);

p_or = (TNr+TPr)/(TNr+FNr+FPr+TPr);
p_er = ((TNr+FPr)*(TNr+FNr)+(FNr+TPr)*(FPr+TPr))/(TNr+FNr+FPr+TPr)^2;
C_kappar = (p_or-p_er)/(1-p_er);

[precr, tprr, fprr, npvr, threshr] = prec_rec(trainOutput, trainTarget);
% -----------
   
save(NNname,'tr','net','cm','cmT','cmV','cmR',...
    'TNt','FNt','FPt','TPt','ACCt',...
    'Pr_Pt','Re_Pt','F1_Pt','Pr_Nt','Re_Nt','F1_Nt','C_kappat',...
    'prect', 'tprt', 'fprt', 'npvt', 'thresht',...
    'TNv','FNv','FPv','TPv','ACCv',...
    'Pr_Pv','Re_Pv','F1_Pv','Pr_Nv','Re_Nv','F1_Nv','C_kappav',...
    'precv', 'tprv', 'fprv', 'npvv', 'threshv',...
    'TN','FN','FP','TP','ACC',...
    'Pr_P','Re_P','F1_P','Pr_N','Re_N','F1_N','C_kappa',...
    'prec', 'tpr', 'fpr', 'npv', 'thresh',...
    'TNr','FNr','FPr','TPr','ACCr',...
    'Pr_Pr','Re_Pr','F1_Pr','Pr_Nr','Re_Nr','F1_Nr','C_kappar',...
    'precr', 'tprr', 'fprr', 'npvr', 'threshr');