%% SKS Data Analyze

%preamble
load('/Volumes/Emily_Data/Seismic_Data_Analysis/sks_data.mat');
addpath '/Volumes/Emily_Data/Seismic_Data_Analysis/mattaup';
addpath '/Volumes/Emily_Data/Seismic_Data_Analysis/ReadMSEEDFast';
addpath '/Volumes/Emily_Data/Seismic_Data_Analysis';
secinday = 24*60*60; 

%parameters
flo = 1/50; %1/longest period
fhi = 1/8; %1/shortest period
dt = 0.02;
sks0 = - 5;
sks1 =  25;
tapertime = 2;


for orid = 19
%for orid = 1:length(eqar)
   
    %filter the data between __ and __ second freq
    [datf] = filt_quick(eqar(orid).datw_sks,flo,fhi,dt);
    
    %plot prewindow, filtered waveform
    figure(3); clf;
    plot(eqar(orid).ttw_sks,datf(:,2),'k');
    hold on;
    xlabel('Time since expected SKS arrivial (sec)');
    ylabel('Radial Component of the Waveform');
    title(['Earthquake number ' num2str(orid) ', magnitude = ' num2str(eqar(orid).eqmag) ', distance = ' num2str(eqar(orid).gcarc)]);
        
    % create __ second window, -__ to +__ relative to sks arrival
    skstime = eqar(orid).skstime;
    check_inwin  = (eqar(orid).ttw_sks >= sks0) & (eqar(orid).ttw_sks < sks1);
    ttw_inwin = eqar(orid).ttw_sks(check_inwin,:);
    datf_inwin = datf(check_inwin,:);
    
    %plot windowed, filtered waveform
    plot(ttw_inwin, datf_inwin(:,2),'g');
    hold on;
     
    %tapering function
    [win] = flat_hanning_taper(ttw_inwin,tapertime);
    win = win(:); %ensures we have a column vector
    datf_taper = datf_inwin .* (win.*ones(1,size(datf_inwin,2)));
    
    %plot windowed, filtered, and tapered wavefrom
    plot(ttw_inwin,datf_taper(:,2),'b');
    hold on;
    legend({'Filtered','Windowed, Filtered','Tapered, Windowed, Filtered'},'Location','northwest');
    
%     chans{1} vertical --> east --> vertical
%     chans{2} north --> --> north --> radial
%     chans{3} east --> vertical --> transverse
    
    % unsplitting function
    datN = eqar(orid).datw_sks_nor;
    datE = eqar(orid).datw_sks_east;
    samprate = 1/eqar(orid).dt;
    plotopt = 1; % option to plot
    %plotopt = 0; % option to not plot

    [phiSC,dTSC, Ru, Tu, Emap, inipolest,phiEV,dTEV,Lmap] = SC_unsplit...
        (datN,datE, samprate, plotopt);
    
%     
%     %plot radial and transverse components against one another
%     figure(5); clf;
%     subplot(2,2,1);
%     plot(datf_taper(:,2));
%     xlim([0 20]);
%     subplot(2,2,2);
%     plot(
%     subplot(2,2,3);
%     plot(datf_taper(:,2),datf_taper(:,3));
%     subplot(2,2,4);
%     plot(Ru, Tu);

end
