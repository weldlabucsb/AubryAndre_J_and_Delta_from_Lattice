%% Quick test of J_Delta_PiecewiseFit

s1s = Outtable3900.('s1');
Js = Outtable3900.('J Numeric');
DeltasWannierIntegral = (2/(beta^2))*Outtable3900.('Delta Numeric');

Des = 0.5*(1064/915)^2*DeltasWannierIntegral;

s2 = 1;

fitJs = zeros(size(s1s));
fitDes = zeros(size(s1s));

gaussJs = zeros(size(s1s));
gaussDes = zeros(size(s1s));

for ii = 1:length(s1s)
    [fitJs(ii), fitDes(ii)] = J_Delta_PiecewiseFit(s1s(ii),s2);
    [gaussJs(ii),gaussDes(ii)] = J_Delta_Gaussian(s1s(ii),s2,1064,915);
end

Jfig = figure();
Jax = axes(Jfig);

plot(Jax,s1s,100*(Js-fitJs)./Js)

title(Jax,'J Residual Error: Piecewise Fit vs Numeric')
ylabel(Jax,'percent error    (100*(numeric - fit)/numeric)')
xlabel(Jax,'s1')


Defig = figure();
Deax = axes(Defig);

plot(Deax,s1s,100*(Des-fitDes)./Des);

title(Deax,'Delta Residual Error: Piecewise Fit vs Numeric')
ylabel(Deax,'percent error    (100*(numeric - fit)/numeric)')
xlabel(Deax,'s1')




JfigGauss = figure();
JaxGauss = axes(JfigGauss);

plot(JaxGauss,s1s,100*(Js-gaussJs)./Js)

title(JaxGauss,'J Residual Error: Gaussian Approx J vs Numeric')
ylabel(JaxGauss,'percent error    (100*(numeric - gaussApprox)/numeric)')
xlabel(JaxGauss,'s1')


DefigGauss = figure();
DeaxGauss = axes(DefigGauss);

plot(DeaxGauss,s1s,100*(Des-gaussDes)./Des);

title(DeaxGauss,'Delta Residual Error: Gaussian Approx Delta vs Numeric')
ylabel(DeaxGauss,'percent error    (100*(numeric - gaussApprox)/numeric)')
xlabel(DeaxGauss,'s1')