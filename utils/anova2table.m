function T = anova2table(ANOVA)
%STATS2TABLE transforms the results structure of the fixed effects coefficients from fixedEffects 
% (e.g., : [beta,betanames,stats] = fixedEffects(lme,'DFMethod','satterthwaite');)
% into a table which could then be saved into an xls file using writetable
% Tamar Regev Jul 25 2023

Name = ANOVA.Term;
FStat = ANOVA.FStat;
DF1 = ANOVA.DF1;
DF2 = ANOVA.DF2;
pValue = ANOVA.pValue;

T = [table(Name),table(FStat), table(DF1), table(DF2),table(pValue)];%



end

