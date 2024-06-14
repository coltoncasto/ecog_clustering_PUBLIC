function T = stats2table(stats)
%STATS2TABLE transforms the results structure of the fixed effects coefficients from fixedEffects 
% (e.g., : [beta,betanames,stats] = fixedEffects(lme,'DFMethod','satterthwaite');)
% into a table which could then be saved into an xls file using writetable
% Tamar Regev Jul 25 2023

Name = stats.Name;
Estimate = stats.Estimate;
SE = stats.SE;
tStat = stats.tStat;
DF = stats.DF;
pValue = stats.pValue;
Lower = stats.Lower;
Upper = stats.Upper;

T = [table(Name),table(Estimate),table(SE),table(tStat), table(DF), table(pValue), table(Lower), table(Upper)];%,tStat,DF,pValue,Lower,Upper];



end

