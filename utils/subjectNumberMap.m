function map=subjectNumberMap

    % ordering of subjects-->participants in paper
    subjects = {'AMC026',...
                'AMC029',...
                'AMC031',...
                'AMC037',...
                'AMC038',...
                'AMC044',...
                'AMC082',...
                'AMC086',...
                'AMC088',...
                'AMC091',...
                'AMC092',...
                'AMC096',...
                'AMC097',...
                'AMC099',...
                'BJH006',...
                'BJH007',...
                'BJH008',...
                'BJH011',...
                'BJH012',...
                'BJH016',...
                'MCJ011',...
                'SLCH002'...
    };
    sessions = num2cell(1:length(subjects));
    map = containers.Map(subjects,sessions);

end