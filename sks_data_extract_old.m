clear all
close all

%% Give top level directory for data card
% This will probably be just /Volumes/Untitled
topdir = '/Volumes/Emily_Data/SR.G1C3'; % no need for final slash
addpath '/Volumes/Emily_Data/SOH_Centaur_Eval/mattaup'
addpath '/Volumes/Emily_Data/SOH_Centaur_Eval/matguts'
addpath '/Volumes/Emily_Data/SOH_Centaur_Eval/IRIS'
addpath '/Volumes/Emily_Data/Seismic_Data_Analysis'
addpath '/Volumes/Emily_Data/Seismic_Data_Analysis/ReadMSEEDFast'

delsac = false; % option to delete SAC files that are made in this dir


%% give station coordinates
slat = 34.69795;
slon = -120.04044;

%% Loop through quakes

load('/Volumes/Emily_Data/SOH_Centaur_Eval/IRIS/events.mat')
    
for orid = 1:length(events_info)
%for orid = 38
    
    %clear eqar;
    elat      = events_info(orid).PreferredLatitude;
    elon      = events_info(orid).PreferredLongitude;
    edep      = events_info(orid).PreferredDepth;
    evtimestr = events_info(orid).PreferredTime;
    eqmag     = events_info(orid).PreferredMagnitudeValue;

    eqar(orid).elat      = elat;
    eqar(orid).elon      = elon;
    eqar(orid).edep      = edep;
    eqar(orid).evtimestr = evtimestr;
    eqar(orid).eqmag     = eqmag;
     

    % no need to change these here
    filtfs = [0.01 0.5];
    viewwind = [-100 2000];% default window to view in seconds from event time
    bigwind = [-500 3500]; % window to subset in seconds from event time

    %% preamble
    addpath('matguts');
    secinday = 24*60*60;

    [data] = make_data_date_struct(topdir);

    %% parse into allfiles
    Ndf = 0; % number of data files
    dfs = {}; % names of data files
    dps = {}; % paths to data files
    serialstarts = []; % serial start times for each file
    for iy = 1:length(data.yrs)
        ystr = num2str(data.yrs(iy),'%04.f');

        for im = 1:data.(['yyyy',ystr])(1).Nm
            mstr = num2str(data.(['yyyy',ystr]).mos(im),'%02.f');
            modat = data.(['yyyy',ystr])(1).(['mm',mstr]);
            Ndf = Ndf + length(modat.datfiles);
            dfs = cat(1,dfs,modat.datfiles);
            dps = cat(1,dps,modat.fpath);
            serialstarts = cat(1,serialstarts,modat.serialstart); %list of start times of all the data files
        end; clear im
    end; clear iy
    

    %% Predict EQ times
    [gcarc,seaz] = distance(slat,slon,elat,elon);
    tpt = tauptime('deg',gcarc,'dep',edep); % measured in seconds, relative time
    
    sks_id = strcmp('SKS',{tpt.phase});
    check1 = any(sks_id);
    if check1 == 0
        fprintf ('%s\n', 'Error: no SKS phase identified');
        continue
    end
        
    evtime = datenum(evtimestr); % measured in days, absolute time
    
    eqar(orid).gcarc  = gcarc;
    eqar(orid).seaz   = seaz;
    eqar(orid).evtime = evtime;
    
    skstime = tpt(sks_id).time;
    sks0 = skstime - 200;
    sks1 = skstime + 200;
    
    % strings to check the dates
    %sks0_str = datestr(sks0);
    %sks1_str = datestr(sks1);
    
    skstime_abs = tpt(sks_id).time/secinday + evtime;
    sks0_abs = skstime_abs - 200/secinday; 
    sks1_abs = skstime_abs + 200/secinday;
 
    % add to eqar structure
    eqar(orid).skstime = skstime;
    eqar(orid).skstime_abs = skstime_abs;
    eqar(orid).sks0 = sks0;
    eqar(orid).sks1 = sks1;
    

    % find data file(s)
    
    if sks1_abs > serialstarts(end) % check if the latest time requested is within the last day
        fprintf ('%s\n', 'Error: earthquake occured outside of data window');
        continue
    end
    
    fid0 = find(serialstarts<sks0_abs,1,'last');
    fid1 = find(serialstarts<sks1_abs,1,'last');
    fid2read = fid0:fid1;
    for idf = 1:length(fid2read)
        msdfile = [dps{fid2read(idf)},'/',dfs{fid2read(idf)}];
        
        % Use ReadMSEEDFast to make a m file in the current directory
        a = ReadMSEEDFast(msdfile);
        chans = unique({a.channel}, 'stable');
        %chans = unique({a.channel});
        tt = [];
        for ii = find(strcmp({a.channel},chans{1}))
          tt = [tt; datenum(a(ii).dateTimeString)+[0:(a(ii).numberOfSamples-1)]'./a(ii).sampleRate./secinday];
        end
        dat = nan(length(tt),length(chans));
        for ic = 1:length(chans)
          dat(:,ic) = double(vertcat(a(strcmp({a.channel},chans{ic})).data));
        end
%         figure(1);clf;
%         for ic = 1:3
%           subplot(3,1,ic); plot(tt_a,dat_a(:,ic))
%         end
    end

    check2 = isempty(fid2read);
    if check2 == 1
        fprintf ('%s\n', 'Error: fid2read was empty');
        continue
    end
    
    eqar(orid).tt = tt;
    dt = round((tt(2) - tt(1))*secinday,3);
    
    
%     figure(1);
%     plot(tt-evtime,dat);
%     plot((tt-evtime)*secinday,filt_quick(dat,1/50,1/8,0.02));
%     xlim ([800 1600]);

   

    % excerpt window of data around event
    inwin_sks = ((tt-evtime)*secinday >= sks0) & ((tt-evtime)*secinday < sks1);
    ttw_sks1 = (tt(inwin_sks,:)-skstime_abs)*secinday;
    datw_sks1 = dat(inwin_sks,:);
%     figure(1); 
%     plot(ttw_sks1,datw_sks1);
%     title ('unfiltered, windowing');
%     figure(2);
%     flo = 1/50; %1/longest period
%     fhi = 1/8; %1/shortest period
%     dt = 0.02;
%     [datf1] = filt_quick(datw_sks1,flo,fhi,dt);
%     plot(ttw_sks1,datf1);
%     title ('filtered, windowing');

    %[datdown] = downsamp(dat,40,10);
  
    %interpolate data onto idealized time vector
    ttw_sks = [sks0:dt:sks1]'-skstime;
    datw_sks = interp1((tt-skstime_abs)*secinday,dat,ttw_sks);
    ttw_sks(end) = []; %delete 40,001 element
    datw_sks(end,:) = []; %delete 40,001 element
%     figure(3);
%     plot(ttw_sks,datw_sks);
%     title ('unfiltered, interpolated');
%     figure(4);
%     [datf] = filt_quick(datw_sks,flo,fhi,dt);
%     plot(ttw_sks,datf);
%     title ('unfiltered, interpolated');
    
%     figure(1);
%     plot(ttw_sks1,datw_sks1,ttw_sks,datw_sks);

%     figure(1);
%     %plot(tt,dat(:,2),ttw_sks1,datw_sks1(:,2),ttw_sks,datw_sks(:,2));
%     plot((tt-skstime_abs)*secinday,dat(:,3)); hold on;
%     plot(ttw_sks1,datw_sks1(:,3)); hold on;
%     plot(ttw_sks,datw_sks(:,3));
%     xlim ([-250 250]);
%     xlabel ('Time since expected SKS arrivial (sec)');
%     ylabel ('Data');
    
    %detrend data
    datw_sks = detrend(datw_sks);
    
    datw_sks = [datw_sks(:,3) datw_sks(:,2) datw_sks(:,1)];
    baz = azimuth(slat,slon,elat,elon);
    [datZRT] = zne2zrt(datw_sks,baz);
       
%   
%     figure(1);
%     clf;
%     %plot((tt-skstime_abs)*secinday,dat(:,2));
%     %hold on;
%     plot(ttw_sks,datZRT(:,2),'k');
%     xlim ([-250 250]);
%     xlabel ('Time since expected SKS arrivial (sec)');
%     ylabel ('Data');
    
    %tapering function
    tapertime = 20;
    [win] = flat_hanning_taper(ttw_sks,tapertime);
    win = win(:); %ensures we have a column vector
    datw_taper = datw_sks .* (win.*ones(1,size(datw_sks,2)));
    
%     figure(1);
%     clf;
%     plot((tt-skstime_abs)*secinday,dat(:,2));
%     hold on;
%     plot(ttw_sks,datw_taper(:,2),'k');
%     xlim ([-250 250]);
%     xlabel ('Time since expected SKS arrivial (sec)');
%     ylabel ('Data');


    % save to eqar structure
    eqar(orid).dt = dt;
    eqar(orid).slat = slat;
    eqar(orid).slon = slon;
    eqar(orid).ttw_sks = ttw_sks;
    eqar(orid).datw_sks = datw_sks;
    eqar(orid).datw_taper = datw_taper;
    
%     save(['/Volumes/Emily_Data/evdata/', events_info(orid).PreferredTime([1:10]),'_',...
%         events_info(orid).PreferredTime([12:13]),'_',events_info(orid).PreferredTime([15:16]),...
%         '_',events_info(orid).PreferredTime([18:19])], 'eqar');
    

end

%delete layers of eqar that are empty/missing information
kill = [];
for i = 1:length(eqar)
    if isempty(eqar(i).dt) 
        kill(i) = 1;
    elseif isempty(eqar(i).sks0)
        kill(i) = 1;
    else
        kill(i) = 0;
    end
end

kill = logical(kill);
eqar(kill) = [];

%save the entire structure
save('sks_data.mat', 'eqar');
