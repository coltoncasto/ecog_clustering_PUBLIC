function T = lme2table(lme)
%LME2TABLE transforms the results structure of the fixed effects coefficients from fitlme into a table which
%could then be saved into an xls file using writetable
% Tamar Regev Dec 19 2020 

Name = lme.Coefficients.Name;
Estimate = lme.Coefficients.Estimate;
SE = lme.Coefficients.SE;
tStat = lme.Coefficients.tStat;
DF = lme.Coefficients.DF;
pValue = lme.Coefficients.pValue;
Lower = lme.Coefficients.Lower;
Upper = lme.Coefficients.Upper;

T = [table(Name),table(Estimate),table(SE),table(tStat), table(DF), table(pValue), table(Lower), table(Upper)];%,tStat,DF,pValue,Lower,Upper];



end

