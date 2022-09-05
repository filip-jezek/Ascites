%% Read data from sample patients

% ProcedureNumber    SubjectID    DateOfProcedure
pacs = readtable('../data/sample paracentesis data for filip 922022.xlsx', ...
    "filetype", 'spreadsheet', 'VariableNamingRule', 'modify', 'Sheet', 'Sheet1');

% ProcedureNumber    SubjectID    DateOfProcedure
procDates = readtable('../data/sample paracentesis data for filip 922022.xlsx', ...
    "filetype", 'spreadsheet', 'VariableNamingRule', 'modify', 'Sheet', 'Sheet2');

% procedure = struct(id, subjectID, date, drainedVolume(x), pressures(x), time(x));
% subject(id, procedure(z))
%%
line = 1;
EOF = false;
subjects = {};
procedures = {};
allProcedures = {};
while ~EOF
    

    sub = str2num(cell2mat(pacs{line, 2}));
    proc = pacs{line, 1};
    nextLineProc = proc;
    vol = [];
    p = [];
    t_abs = [];
    % assuming it is sorted by proc
    while nextLineProc == proc

        vol = [vol pacs{line, 3}];
        p = [p pacs{line, 4}];
        try
            t = datetime(pacs{line, 5}, 'InputFormat', 'hh:mma');
        catch
            t = NaT;
        end
        t_abs = [t_abs t];

        if line + 1 > size(pacs, 1)
            EOF = true;
            break;
        else
            line = line + 1;
        end
        nextLineProc = pacs{line, 1};
        nextSub = str2num(cell2mat(pacs{line, 2}));
    end
   
    % find procedure time in second spreadsheet
    procDate = NaT;
    for pd = 1:size(procDates, 1)
        if str2num(cell2mat(procDates{pd, 1})) == proc
            procDate = datetime(procDates{pd, 3});break;
        end
    end
    
    %     times = minutes(t_abs - t_abs(1));
    times = timeofday(t_abs) + procDate;
    procedure = struct('Id', proc, 'SubjectId', sub, 'Volumes', vol, ...
        'Pressures', p, 'Times', times, 'Date', procDate);
    
    
    % look if we do not already have this subject
    theresheis = false;
    for s = 1:length(subjects)
        if subjects{s}.Id == sub 
            theresheis = true;
            break;
        end
    end
    if theresheis
        subjects{s}.Procedures{end + 1} = procedure;
    else
        procedures{1} = {procedure};
        subject = struct('Id', sub, 'Procedures', procedures);
        subjects{end + 1} = subject;
    end
    
    allProcedures{end + 1} = procedure; 
end


%% All subjects separately
figure(1);clf;

for s = 1:15
  
subplot(3,5,s);hold on;
prc = subjects{s}.Procedures;
    for p = 1:length(prc)
        plot(prc{p}.Volumes, prc{p}.Pressures);
    end
    title(['Subject ' num2str(s) ' with ' num2str(length(prc)) ' prc']);
    xlim([0, 10]);
    ylim([0, 20]);
%     set(gca, 'XDir','reverse');
    if s == 1 || s == 6 || s == 11
        ylabel('Pressure (mmHg)')
    end
    if s > 10
        xlabel('Drained volume (L)');
    end
end

%% All subjects at once
figure(2);clf; hold;
for s = 1:15
  
prc = subjects{s}.Procedures;
    for p = 1:length(prc)
        plot(prc{p}.Volumes, prc{p}.Pressures);
    end
end
xlim([0, 10]);
ylim([0, 20]);
ylabel('Pressure (mmHg)')
xlabel('Drained volume (L)');

%% All subjects at once - max volumes and pressures only
figure(2);clf; hold on;
subnum = length(subjects);
cols = turbo(subnum);
markers = repmat(['s'; 'o'; 'v'], [ceil(subnum/3), 1])
for s = 1:subnum
    prc = subjects{s}.Procedures;
    for p = 1:length(prc)
       plo(s) = plot(max(prc{p}.Volumes), max(prc{p}.Pressures), markers(s), 'Color', cols(s, :), 'MarkerFaceColor', cols(s, :));
    end
    l{s} = ['Subject ' num2str(s)];
end
xlim([0, 15]);
ylim([0, 30]);
ylabel('Pressure (mmHg)')
xlabel('Drained volume (L)');
legend(plo, l);
title('Maximal drained volume');

%% Readmission dates
for sub = 1:length(subjects)
    if length(subjects{sub}.Procedures) < 2
        continue;
    end
    for lebn
    pv = max(subjects{sub}.Procedures.vol); % drained volume
    pp = max(subjects{sub}.Procedures.p); % starting pressure
    prp = min(subjects{sub}.Procedures.p); % resting (ending) pressure
    pst = subjects{sub}.Procedures.Times(0); % start time
    pet = subjects{sub}.Procedures.Times(end); % end time
    
end