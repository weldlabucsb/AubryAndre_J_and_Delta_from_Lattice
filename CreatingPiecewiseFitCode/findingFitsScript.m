%% Script for Fitting the Numerical values J and Delta with 1064 and 915
%All this is following formulas in Michele Modugno 2009 New J. Phys. 11, 033023

la1 = 1064;
la2 = 915;

beta = la1/la2;

% Outtable3900 is the output from compare_J_Delta_Methods.m. 
%  Named Outtable3900 because there were 3901 points for this fit
s1s = Outtable3900.('s1');
Js = Outtable3900.('J Numeric');
DeltasWannierIntegral = (2/(beta^2))*Outtable3900.('Delta Numeric');


% Define fittype
JFitType = fittype('A*s1^alp*exp(-B*s1^beta)','independent',{'s1'},...
    'coefficients',{'A','alp','B','beta'},'dependent','J');
DeltaFitType = fittype('A*exp(-C*s1^(-gamma))','independent',{'s1'},...
    'coefficients',{'A','C','gamma'});

% Define Fit options
JFitOpts = fitoptions(JFitType);
JFitOpts.StartPoint = [1.43,0.98,2.07,0.5];

DeltaFitOpts = fitoptions(DeltaFitType);
DeltaFitOpts.StartPoint = [1,beta^2,0.5];

% Define Fits
JFit = fit(s1s,Js,JFitType,JFitOpts);
DeltaFit = fit(s1s,DeltasWannierIntegral,DeltaFitType,DeltaFitOpts);


% Check Fit goodness

f_J = figure();
ax_J = axes();
hold(ax_J,'on');

plot(ax_J,s1s,Js)
plot(ax_J,s1s,JFit(s1s))

title(ax_J,'J vs s1')



f_JResidue = figure();
ax_JResidue = axes();
hold(ax_JResidue,'on');

plot(ax_JResidue,s1s,100*(Js-JFit(s1s))./Js)

title(ax_JResidue,'J Residuals')
ylabel(ax_JResidue,'percent error')
xlabel(ax_JResidue,'s1')


%% J Fits
% 17 to 40Er range
s1sOver17 = s1s(s1s>=17);
JsOver17 = Js(s1s>=17);

JFitOver17 = fit(s1sOver17,JsOver17,JFitType,JFitOpts);

f_JResidueOver17 = figure();
ax_JResidueOver17 = axes();
hold(ax_JResidueOver17,'on');

plot(ax_JResidueOver17,s1sOver17,100*(JsOver17-JFitOver17(s1sOver17))./JsOver17)

title(ax_JResidueOver17,'J Residuals s1 >= 17')
ylabel(ax_JResidueOver17,'percent error')
xlabel(ax_JResidueOver17,'s1')

% 12 to 17Er range

s1s12to17 = s1s((s1s<=17)&(s1s>=12));
Js12to17 = Js((s1s<=17)&(s1s>=12));

JFit12to17 = fit(s1s12to17,Js12to17,JFitType,JFitOpts);

f_JResidue12to17 = figure();
ax_JResidue12to17 = axes();
hold(ax_JResidue12to17,'on');

plot(ax_JResidue12to17,s1s12to17,100*(Js12to17-JFit12to17(s1s12to17))./Js12to17)

title(ax_JResidue12to17,'J Residuals s1 12 to 17Er')
ylabel(ax_JResidue12to17,'percent error')
xlabel(ax_JResidue12to17,'s1')


% 8 to 12Er range

s1s8to12 = s1s((s1s<=12)&(s1s>=8));
Js8to12 = Js((s1s<=12)&(s1s>=8));

JFit8to12 = fit(s1s8to12,Js8to12,JFitType,JFitOpts);

f_JResidue8to12 = figure();
ax_JResidue8to12 = axes();
hold(ax_JResidue8to12,'on');

plot(ax_JResidue8to12,s1s8to12,100*(Js8to12-JFit8to12(s1s8to12))./Js8to12)

title(ax_JResidue8to12,'J Residuals s1 8 to 12Er')
ylabel(ax_JResidue8to12,'percent error')
xlabel(ax_JResidue8to12,'s1')


% 4 to 8Er range  POLY FIT

s1s4to8 = s1s((s1s<=8)&(s1s>=4));
Js4to8 = Js((s1s<=8)&(s1s>=4));

JFit4to8 = fit(s1s4to8,Js4to8,'poly3');

f_JResidue4to8 = figure();
ax_JResidue4to8 = axes();
hold(ax_JResidue4to8,'on');

plot(ax_JResidue4to8,s1s4to8,100*(Js4to8-JFit4to8(s1s4to8))./Js4to8)

title(ax_JResidue4to8,'J Residuals s1 4 to 8Er')
ylabel(ax_JResidue4to8,'percent error')
xlabel(ax_JResidue4to8,'s1')


% 1 to 4Er range  POLY FIT

s1s1to4 = s1s((s1s<=4)&(s1s>=1));
Js1to4 = Js((s1s<=4)&(s1s>=1));

JFit1to4 = fit(s1s1to4,Js1to4,'poly4');

f_JResidue1to4 = figure();
ax_JResidue1to4 = axes();
hold(ax_JResidue1to4,'on');

plot(ax_JResidue1to4,s1s1to4,100*(Js1to4-JFit1to4(s1s1to4))./Js1to4)

title(ax_JResidue1to4,'J Residuals s1 1 to 4Er')
ylabel(ax_JResidue1to4,'percent error')
xlabel(ax_JResidue1to4,'s1')


%% Delta Wannier Integral Fits

figDelta=figure();
axDelta = axes(figDelta);
plot(axDelta,s1s,DeltasWannierIntegral);
title(axDelta,'Delta Wannier Integral as a function of s1')
xlabel(axDelta,'s1')
ylabel(axDelta,'Delta Integral (aka Delta / (s2*beta^2/2))')


% 17 to 40Er range
s1sOver17 = s1s(s1s>=17);
DeltaIntsOver17 = DeltasWannierIntegral(s1s>=17);

DeltaIntsFitOver17 = fit(s1sOver17,DeltaIntsOver17,DeltaFitType,DeltaFitOpts);

f_DIntResidueOver17 = figure();
ax_DIntResidueOver17 = axes();
hold(ax_DIntResidueOver17,'on');

plot(ax_DIntResidueOver17,s1sOver17,100*(DeltaIntsOver17-DeltaIntsFitOver17(s1sOver17))./DeltaIntsOver17)

title(ax_DIntResidueOver17,'Delta Wannier Integral Residuals s1 >= 17')
ylabel(ax_DIntResidueOver17,'percent error')
xlabel(ax_DIntResidueOver17,'s1')


% 8 to 17Er range
s1s8to17 = s1s((s1s<=17)&(s1s>=8));
DeltaInts8to17 = DeltasWannierIntegral((s1s<=17)&(s1s>=8));

DeltaIntsFit8to17 = fit(s1s8to17,DeltaInts8to17,DeltaFitType,DeltaFitOpts);

f_DIntResidue8to17 = figure();
ax_DIntResidue8to17 = axes(f_DIntResidue8to17);
hold(ax_DIntResidue8to17,'on');

plot(ax_DIntResidue8to17,s1s8to17,100*(DeltaInts8to17-DeltaIntsFit8to17(s1s8to17))./DeltaInts8to17)

title(ax_DIntResidue8to17,'Delta Wannier Integral Residuals s1 8 to 17')
ylabel(ax_DIntResidue8to17,'percent error')
xlabel(ax_DIntResidue8to17,'s1')


% 4 to 8Er range
s1s4to8 = s1s((s1s<=8)&(s1s>=4));
DeltaInts4to8 = DeltasWannierIntegral((s1s<=8)&(s1s>=4));

DeltaIntsFit4to8 = fit(s1s4to8,DeltaInts4to8,DeltaFitType,DeltaFitOpts);

f_DIntResidue4to8 = figure();
ax_DIntResidue4to8 = axes(f_DIntResidue4to8);
hold(ax_DIntResidue4to8,'on');

plot(ax_DIntResidue4to8,s1s4to8,100*(DeltaInts4to8-DeltaIntsFit4to8(s1s4to8))./DeltaInts4to8)

title(ax_DIntResidue4to8,'Delta Wannier Integral Residuals s1 4 to 8')
ylabel(ax_DIntResidue4to8,'percent error')
xlabel(ax_DIntResidue4to8,'s1')

% 1 to 4Er range
s1s1to4 = s1s((s1s<=4)&(s1s>=1));
DeltaInts1to4 = DeltasWannierIntegral((s1s<=4)&(s1s>=1));

DeltaIntsFit1to4 = fit(s1s1to4,DeltaInts1to4,'poly4');

f_DIntResidue1to4 = figure();
ax_DIntResidue1to4 = axes(f_DIntResidue1to4);
hold(ax_DIntResidue1to4,'on');

plot(ax_DIntResidue1to4,s1s1to4,100*(DeltaInts1to4-DeltaIntsFit1to4(s1s1to4))./DeltaInts1to4)

title(ax_DIntResidue1to4,'Delta Wannier Integral Residuals s1 1 to 4')
ylabel(ax_DIntResidue1to4,'percent error')
xlabel(ax_DIntResidue1to4,'s1')
