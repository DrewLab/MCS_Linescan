function [mv_mpP]=CalculateVelocityLinescanTIFFMCS(sampling_freq,N_analog_channels)
% opens a tiff file exported from mpview 2.0, gets xdimension to scan over from user input, makes window size
% calculates velocity from linescans using the Radon transform
% sampling_freq - (in Hz), Should be >30 to avoid aliasing mouse heart rate
% N_analog_channels - # of analog channels to read into analog data, if 0,
% don't read in any files

% structures are mv_mpP, with mpP substructs
% pjd 5/2017
% Note lots of legacy code to deal with non-interleaved line-scans
% requires Chronux
my_crunched_folder='/Users/patrick/Documents/data/crunched/';% change this to your own crunched directory
matfilename=[];%contains arbitraryscanpath data
line_boundaries=[];
[filename,pathname]=uigetfile('*.TIF','pick a linescan file')%loads file



if N_analog_channels==1%only works with 1 channel for now
    [analog_data_filename,analog_data_pathname,]=uigetfile('*.txt','pick a text file with analog data')%loads file
    cd(analog_data_pathname);
    analog_data=load(analog_data_filename,'ascii');
end
cd(pathname)% go to the data file
%display the file
[xstart, xstop,~,~,vessel_ID,depths,scan_type]=Display_Frames_MultiVessel_TIFF(filename,matfilename);
nvessels=length(xstart);%how many vessels?
startframe=1;%start at the beginning  %input('start frame: ');
endframe=length(imfinfo(filename))%start at the end input('end frame: ');
the_decimate=1;%input('decimation factor: '); only for oversampled data
user_microns_per_pixel=input('microns per pixel @ 1x');
the_objective=input('objective: [1]-40x Olympus: [2]-20x, 0.5NA: [3]-20x, 0.95NA: [4]-10x: [5]-4x: [6]-other:');
animal=input('animal?','s');

max_v=10000;% maximum velocity 
windowsize=1/(sampling_freq);
angle_span=15;%number of degrees (+/-) around previous angle to look
angle_resolution=.1;%how accurate to determine the angle
the_scanmirrors=2;%we only use 6215 input('which mirrors? 1)6210(fast) 2)6215(slow)');
channel=1;

for v=1:nvessels
    thisfile=mfilename('fullpath');
    mpP.code=GetFunctionCode(thisfile);
    if (scan_type(v)==1)%get velocity
        [mpP,function_ID]=GetVelocityRadon_MCS_L4_cc_TIFF(filename,...
            [startframe endframe],windowsize,angle_span, angle_resolution,[xstart(v) xstop(v)],...
            the_decimate,max_v,user_microns_per_pixel,channel,the_scanmirrors);%
        
        try
            mpP.radon_code=GetFunctionCode(function_ID);
            
        catch
            disp('code read fail!')
        end
    end
    
    
    if the_scanmirrors==1% legacy
        mpP.the_scan_mirrors=6210;
    elseif the_scanmirrors==2% legacy
        mpP.the_scan_mirrors=6215;
    end
    try
        mpP.ECoG.raw=mpP.Ch3;
    end
    mpP.animal=animal;
    
    try
        mpP.vessel=vessel_ID{v};
    end
    mpP.clock=clock;
    mpP.the_objective=the_objective;
    mpP.stimtime=NaN;
    mpP.depth=depths(v);
    try
        mpP.scan_type=scan_type(v);
    end
    if v>1
        [mv_mpP,mpP]=MatchStructs(mv_mpP,mpP);
        mv_mpP=orderfields(mv_mpP);
        mpP=orderfields(mpP);
    end
    mv_mpP(v)=mpP;%put the substruct in the structure
    try
        mv_mpP(1).analog_data=analog_data;%put in the analog data
        mv_mpP(1).analog_data_filename=analog_data_filename;
    catch
        mv_mpP(1).analog_data=[];%put in the analog data
    end
end
ntrials=1;%input
the_date=date;
[mv_mpP]=linescan_xcov_velocity_TIFF(mv_mpP);%calculate the velocity using the cross correlation method of Kim et al, 2012

save_filename=[animal '_' vessel_ID{1} '_mv_mpP_' filename(1:end-4) '_' the_date];

try
    save(save_filename,'mv_mpP')%save processed file to current directory
    save([my_crunched_folder save_filename],'mv_mpP')%save to crunched directory
catch
    disp('crunched folder save failed!')
end
try
    plot_mv_mpP_TIFF_MCS(mv_mpP,109,'all');
catch
    disp('plot failed')
end
end

%% subfunctions
function [the_code]=GetFunctionCode(filename)
%this function the file and returns it as a string
a = fopen([ filename '.m']);
the_code_numbers=fread(a);
the_code = char(the_code_numbers);
end

function [s1,s2]=MatchStructs(s1,s2)
fields1=fieldnames(s1);
fields2=fieldnames(s2);
mpP=s2;
for n=1:length(fields1)
    if (strcmp(fields1{n},fields2)==0)
        for k=1:length(s2)
            s2(k).(fields1{n})=[];%=setfield(s2(k),fields1{n},[]);%put an empty field in
        end
    end
end

for n=1:length(fields2)
    if (strcmp(fields2{n},fields1)==0)
        for k=1:length(s1)
            s1(k).(fields2{n})=[];%s1(k)=setfield(s1(k),fields2{n},[]);%put an empty field in
        end
    end
end
end

function [xstart, xstop,frames_hold,nframes,vessel_ID,depths,scan_type]=Display_Frames_MultiVessel_TIFF(fname,mfname)
%this function displays the first frames of the TIFF linescan file and gets the
%user input to determine the x range to use for the radon transform
%Read header and take further action based on header information
% gives the user the opportunity to record vessel ID and depths
[mpP.Header]=ReadMCSTiffInfo(fname);

nframes=mpP.Header.Frame_Count;
mpP.xsize = str2double(mpP.Header.Frame_Width);
mpP.ysize = str2double(mpP.Header.Frame_Height);
frame_height=mpP.ysize;

frames_hold=LoadTiffConcatenate(fname,[1 5]);%
figure(2)
colormap gray
%plot the data
hold off
imagesc(double(frames_hold))
axis image
hold on
%
axis off
title(fname)
%plot the linescan trajectory
if length(mfname)>0 %plot the linescan trajectory
    try
        load(mfname)
        line_boundaries=find(abs(diff(scanData.pathObjSubNum)));
        scan_velocity=scanData.scanVelocity;
        for ss=1:length(line_boundaries)
            plot(line_boundaries(ss)*[1 1],[0 1000],'r')
        end
        subplot(212)
        hold off
        imagesc(scanData.im)
        axis image
        hold on
        pathImCoords(:,1) = scanData.path(:,1) * (size(scanData.im,2)-1)/(abs(diff(scanData.axisLimCol))) + 1 - (size(scanData.im,2)-1)/(abs(diff(scanData.axisLimCol)))*min(scanData.axisLimCol);
        pathImCoords(:,2) = scanData.path(:,2) * (size(scanData.im,1)-1)/(abs(diff(scanData.axisLimRow))) + 1 - (size(scanData.im,1)-1)/(abs(diff(scanData.axisLimRow)))*min(scanData.axisLimRow);
        plot(pathImCoords(:,1),pathImCoords(:,2),'r')
        
    end
end
nvessels=input('How many vessels?');
scan_type=-ones(nvessels,1);
vessel_ID=cell(nvessels,1);
depths=-ones(nvessels,1);

for vessel=1:nvessels
    disp(['select left and right bounds for Radon analysis of  vessel #' num2str(vessel)])
    [x,y]=ginput(2);
    xstart(vessel)=round(min(x))
    xstop(vessel)=round(max(x))
    depths(vessel)=input('depth?')
    vessel_ID{vessel}=input('vessel ID?','s')
    scan_type(vessel)=input('Scan type? 0=diameter 1=velocity')
end
end

function [mpP, fixed_v_noHR]=HeartRateRemove(mpP)
%removes high freqency (>8hz) noise with 2nd order butterworth filter
f_high=5;
fs=mpP.Blood_flow.the_Fs;
[hr_B,hr_A]=butter(2,2*[f_high]/fs,'low');
for k=1:min(size(mpP.Blood_flow.fixed_v))
    fixed_v_noHR(k,:)=filter(hr_B,hr_A,mpP.Blood_flow.fixed_v(k,:));%filter the velocity signal
end
mpP.Blood_flow.fixed_v_noHR=fixed_v_noHR;
end

function [mpP,function_ID]=GetVelocityRadon_MCS_L4_cc_TIFF(fname,theframes,window_time_size,angle_span, angle_resolution,the_x,the_decimate, max_v, user_microns_per_pixel,channel,isinterleaved);
%Modified version of GetVelocityRadon_MPScan_L4_cc to work with tiff files from the Sutter MCS program.
%PJD 5/2017

%filename,matfilename,...
%            [startframe endframe],windowsize,angle_span, angle_resolution,[xst%art(v) xstop(v)],...
%            the_decimate,max_v,user_microns_per_pixel,channel,the_scanmirrors)%funtion for loading *.MPD linescan files and performing the Radon transform on
%them in order to get the RBC velocity
%
%becase MCS gets line scans in both directions, the frame is split into even
%and odd lines and velocity is extravted independently from each one and then averaged
%INPUTS:
%fname-filename
%matfilename-matalab file containg scan path etc
%theframes - beginning and ending frames
%windowtimesize - size (in sec) of the advancemt of the Radon window
%anglespan - how many degrees (+/-) to look in the adaptive search
%angle_resolution - fine angle resolution
%the_x - [left  right] bounds of which to extract from the frames in the X
%dimension
%the_decimate - factor by which to decimate the line scan in the X dimension
%max_v - time points above this absolute speed are clipped out and interpolated
%between the last and next valid time points
%the_objective - a number designating which objective is used isinterleaved
%user_microns_per_pixel-at 1x for current objective+microscope setup

%OUTPUTS
%thetas - angles of forward and back linescans
%V_hold - velocities of forward and back linescans
%v_out - 0.5*(v_forward+v_back)
%mpP -  a structure with the relevant data
%mpP.Header substruct
%mpP.Ch1 first two frames from channel 1
%mpP.Blood_flow structure with velocity info etc.


%Open stream
function_ID = mfilename('fullpath');

rawinfo=imfinfo(fname);%get the image info
[Header]=ReadMCSTiffInfo(fname);
mpP.Header=Header;

%Read header and take further action based on header information
mpP.num_frames = theframes(end)-theframes(1);%str2double(mpContent.Header.Frame_Count);
mpP.xsize = str2double(mpP.Header.Frame_Width);
mpP.ysize = str2double(mpP.Header.Frame_Height);
mpP.xstart=the_x(1);
mpP.xstop=the_x(end);
frame_height=mpP.ysize;



if strcmp(mpP.Header.Scan_Mode,'Line Scan')
    time_per_line=(mpP.Header.Frame_Width)*(5/4)*(mpP.Header.Pixel_Clock)*(.05*1e-6);
    isinterleaved=1;
end

% if strcmp(mpP.Header.Scan_Mode,'Arbitrary Scan')%
%     time_per_line=(mpP.Header.Frame_Width)*(mpP.Header.Pixel_Clock)*(.05*1e-6);
%     isinterleaved=0; %arb scans are not interleace
% end
mpP.Header.LineRate=1/time_per_line;
Fs_Blood_Flow=mpP.Header.LineRate;
windowsizea=round(Fs_Blood_Flow*window_time_size);%*.25 %for interlevaed lines

%adapt windowsize to interleaved/not interleaved
if isinterleaved==1
    windowsize=8*round(windowsizea/8);
else
    windowsize=4*round(windowsizea/4);
end
mpP.Blood_flow.windowsize=windowsize;


%%%%basic stuff
nframes=theframes(2)-theframes(1)+1;
stepsize=.25*windowsize;
nlines=nframes*(mpP.Header.Frame_Height);
npoints=max(the_x)-min(the_x)+1;%str2num(mpP.Header.Frame_Width);

nsteps=floor(nlines/stepsize)-3;
angles=(1:180);%parameters for initial search
angles_baseline=1:15:180;
angles_adaptive=[-angle_span:1:angle_span]; %angle range to look over
angles_fine=-1.75:angle_resolution:1.75; % fine grained search
n_adaptive_angles=length(angles_adaptive);
theta_low=5;%boundaries for for adaptivce theta
theta_high=n_adaptive_angles-5;
spread_matrix=zeros(1,length(angles));
spread_matrix_adaptive=zeros(2,nsteps,length(angles_adaptive));
baseline_var=zeros(2,nsteps,length(angles_baseline));
spread_matrix_fine=zeros(2,nsteps,length(angles_fine));
thetas=zeros(2,nsteps);%the angle of the lines
data_var=zeros(nsteps,1);
data_hold=zeros(windowsize,npoints);

%data caching
framerate=round(mpP.Header.LineRate/mpP.ysize);
nframes_to_cache=round(mpP.Header.LineRate/mpP.ysize);
data_temp=LoadTiffConcatenate(fname,theframes);%load froma tiff file
data_cache=data_temp(:,the_x(1):the_decimate:the_x(2)-1);
mpP.Blood_flow.mean_BG=(mean(data_cache));
cached_lines=1:(nframes_to_cache*mpP.ysize);
use_lines=1:windowsize;

%
if isinterleaved==1
    data_hold=zeros(2,windowsize/2,length([the_x(1):the_decimate:the_x(2)-1]));
    data_hold(1,:,:)=double(data_cache(1:2:windowsize,1:the_decimate:end));
    data_hold(2,:,:)=double(data_cache(2:2:windowsize,end:-the_decimate:1));
else %not used
    data_hold=zeros(1,windowsize,length([the_x(1):the_decimate:the_x(2)-1]));
    data_hold(1,:,:)=double(data_cache(1:windowsize,1:the_decimate:end));
end
line_counter=4*stepsize;%this is the end of the window

data_cachemean_forward=mean(data_cache);
if isinterleaved==1
    data_cachemean_backward=(fliplr(data_cachemean_forward));
end

%%
if isinterleaved==1
    for k=1:nsteps
        frame_line_marker=line_counter;%find the new place in the frame
        data_hold(1,:,:)=double(data_cache(frame_line_marker-4*stepsize+1:2:frame_line_marker-1,1:the_decimate:end));
        data_hold(2,:,:)=double(data_cache(frame_line_marker+2-4*stepsize:2:frame_line_marker,end:-the_decimate:1));
        use_lines=use_lines+stepsize;
        %take out the local mean of the window
        the_t(k)=1+(k-1)*stepsize+windowsize/2;
        npoints_decimated=size(data_hold,3);
        for n=1:npoints_decimated
            data_hold_ms(1,:,n)=data_hold(1,:,n)-data_cachemean_forward(n);
            data_hold_ms(2,:,n)=data_hold(2,:,n)-data_cachemean_backward(n);
        end
        data_hold_ms=data_hold_ms-mean(data_hold_ms(:));
        if k==1
            radon_hold(1,:,:)=(radon(squeeze(data_hold_ms(1,:,:)),angles));
            radon_hold(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)),angles);
            spread_matrix(1,:)=var(squeeze(radon_hold(1,:,:)));
            spread_matrix(2,:)=var(squeeze(radon_hold(2,:,:)));
            [m(1,:) the_theta(1,:)]=max(spread_matrix(1,k));
            [m(2,:) the_theta(2,:)]=max(spread_matrix(2,k));
            thetas(1,k)=angles(the_theta(1,:));
            thetas(2,k)=angles(the_theta(2,:));
            radon_hold_adaptive(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)),-thetas(1,k)+angles_adaptive);
            radon_hold_adaptive(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)),-thetas(2,k)+angles_adaptive);
            spread_matrix_adaptive(1,k,1)=max(var(squeeze(radon_hold_adaptive(1,:,:))));
            spread_matrix_adaptive(2,k,1)=max(var(squeeze(radon_hold_adaptive(2,:,:))));
            radon_hold_fine(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)-mean(data_hold_ms(:))),thetas(1,k)+angles_fine);
            radon_hold_fine(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)-mean(data_hold_ms(:))),thetas(2,k)+angles_fine);
            spread_matrix_fine(1,k,:)=squeeze(var(squeeze(radon_hold_fine(1,:,:))));
            spread_matrix_fine(2,k,:)=squeeze(var(squeeze(radon_hold_fine(2,:,:))));
            [m the_theta(1,:)]=max(squeeze(spread_matrix_fine(1,k,:)));
            [m the_theta(2,:)]=max(squeeze(spread_matrix_fine(2,k,:)));
            thetas(1,k)=thetas(1,k)+angles_fine(the_theta(1,:));
            thetas(2,k)=thetas(2,k)+angles_fine(the_theta(2,:));
            baseline_var(1,k,:)=var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));
            baseline_var(2,k,:)=var(radon(squeeze(data_hold_ms(2,:,:)),angles_baseline));
        else
            thetas(1,k)=thetas(1,k-1);
            thetas(2,k)=thetas(2,k-1);
            radon_hold_adaptive(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)),thetas(1,k-1)+angles_adaptive);
            radon_hold_adaptive(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)),thetas(2,k-1)+angles_adaptive);
            spread_matrix_adaptive(1,k,:)=var(squeeze(radon_hold_adaptive(1,:,:)));
            spread_matrix_adaptive(2,k,:)=var(squeeze(radon_hold_adaptive(2,:,:)));
            [m the_theta(1,:)]=max(spread_matrix_adaptive(1,k,:));
            [m the_theta(2,:)]=max(spread_matrix_adaptive(2,k,:));
            baseline_var(1,k,:)=var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));%figure out the average radon of the section for calculating S/N
            baseline_var(2,k,:)=var(radon(squeeze(data_hold_ms(2,:,:)),angles_baseline));%figure out the average radon of the section for calculating S/N
            %if the peak of the variance is at an extreme, redo for angles centered
            %around the peak
            
            while ((the_theta(1,:)<=theta_low)||(the_theta(1,:)>=theta_high))
                thetas(1,k)=thetas(1,k)+angles_adaptive(the_theta(1,:));
                radon_hold_adaptive(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)), thetas(1,k)+angles_adaptive);
                spread_matrix_adaptive(1,k,:)=var(squeeze(radon_hold_adaptive(1,:,:)));
                [m the_theta(1,:)]=max(squeeze(spread_matrix_adaptive(1,k,:)));
            end
            while ((the_theta(2,:)<=theta_low)||(the_theta(2,:)>=theta_high))
                thetas(2,k)=thetas(2,k)+angles_adaptive(the_theta(2,:));
                radon_hold_adaptive(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)), thetas(2,k)+angles_adaptive);
                spread_matrix_adaptive(2,k,:)=var(squeeze(radon_hold_adaptive(2,:,:)));
                [m the_theta(2,:)]=max(squeeze(spread_matrix_adaptive(2,k,:)));
            end
            
            baseline_var(1,k,:)=var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));%
            baseline_var(2,k,:)=var(radon(squeeze(data_hold_ms(2,:,:)),angles_baseline));%
            thetas(1,k)=thetas(1,k)+angles_adaptive(the_theta(1,:));
            radon_hold_fine(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)),thetas(1,k)+angles_fine);
            spread_matrix_fine(1,k,:)=var(squeeze(radon_hold_fine(1,:,:)));
            [m the_theta(1,:)]=max(spread_matrix_fine(1,k,:));
            thetas(1,k)=thetas(1,k)+angles_fine(the_theta(1,:));
            thetas(2,k)=thetas(2,k)+angles_adaptive(the_theta(2,:));
            radon_hold_fine(2,:,:)=radon(squeeze(data_hold_ms(2,:,:)),thetas(2,k)+angles_fine);
            spread_matrix_fine(2,k,:)=var(squeeze(radon_hold_fine(2,:,:)));
            [m the_theta(2,:)]=max(spread_matrix_fine(2,k,:));
            thetas(2,k)=thetas(2,k)+angles_fine(the_theta(2,:));
        end
        line_counter=line_counter+stepsize;
    end
end
%% not used
if isinterleaved==0 %non interleved data
    % not yet working with MCS software
    %     for k=1:nsteps
    %         line_counter=line_counter+stepsize;
    %         frame_line_marker=find(cached_lines==line_counter+1);%find the new place in the frame
    %         %%%check to see if the data hold will overrun the cache
    %
    %         data_hold(1,:,:)=double(data_cache(frame_line_marker:frame_line_marker+4*stepsize-1,1:the_decimate:end-1));
    %
    %         use_lines=use_lines+stepsize;
    %         the_t(k)=1+(k-1)*stepsize+windowsize/2;
    %         npoints_decimated=size(data_hold,3);
    %         for n=1:npoints_decimated
    %             data_hold_ms(1,:,n)=data_hold(1,:,n)-data_cachemean_forward(n);
    %         end
    %         data_hold_ms=data_hold_ms-mean(data_hold_ms(:));
    %         if k==1
    %             radon_hold(1,:,:)=(radon(squeeze(data_hold_ms(1,:,:)),angles));
    %             spread_matrix(1,:)=var(squeeze(radon_hold(1,:,:)));
    %             [m(1,:) the_theta(1,:)]=max(spread_matrix(1,k));
    %             thetas(1,k)=angles(the_theta(1,:));
    %             radon_hold_adaptive(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)),-thetas(1,k)+angles_adaptive);
    %             spread_matrix_adaptive(1,k,1)=max(var(squeeze(radon_hold_adaptive(1,:,:))));
    %             radon_hold_fine(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)-mean(data_hold_ms(:))),thetas(1,k)+angles_fine);
    %             spread_matrix_fine(1,k,:)=squeeze(var(squeeze(radon_hold_fine(1,:,:))));
    %
    %             [m the_theta(1,:)]=max(squeeze(spread_matrix_fine(1,k,:)));
    %             thetas(1,k)=thetas(1,k)+angles_fine(the_theta(1,:));
    %             baseline_var(1,k,:)=var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));
    %         else
    %             thetas(1,k)=thetas(1,k-1);
    %             radon_hold_adaptive(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)),thetas(1,k-1)+angles_adaptive);
    %             spread_matrix_adaptive(1,k,:)=var(squeeze(radon_hold_adaptive(1,:,:)));
    %             [m the_theta(1,:)]=max(spread_matrix_adaptive(1,k,:));
    %             baseline_var(1,k,:)=var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));%figure out the average radon of the section for calculating S/N
    %
    %             while ((the_theta(1,:)<=theta_low)||(the_theta(1,:)>=theta_high))
    %                 thetas(1,k)=thetas(1,k)+angles_adaptive(the_theta(1,:));
    %                 radon_hold_adaptive(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)), thetas(1,k)+angles_adaptive);
    %                 spread_matrix_adaptive(1,k,:)=var(squeeze(radon_hold_adaptive(1,:,:)));
    %                 [m the_theta(1,:)]=max(squeeze(spread_matrix_adaptive(1,k,:)));
    %             end
    %             baseline_var(1,k,:)=var(radon(squeeze(data_hold_ms(1,:,:)),angles_baseline));%
    %             thetas(1,k)=thetas(1,k)+angles_adaptive(the_theta(1,:));
    %             radon_hold_fine(1,:,:)=radon(squeeze(data_hold_ms(1,:,:)),thetas(1,k)+angles_fine);
    %             spread_matrix_fine(1,k,:)=var(squeeze(radon_hold_fine(1,:,:)));
    %             [m the_theta(1,:)]=max(spread_matrix_fine(1,k,:));
    %             thetas(1,k)=thetas(1,k)+angles_fine(the_theta(1,:));
    %
    %         end
    %     end
end
%end of radon loop

%%

if((length(find(thetas<0))>=1)&&(length(find(thetas>0))>=1))
    disp('sign flip!') 
end

size(spread_matrix_adaptive);
mpP.Blood_flow.thetas=thetas;

mpP.Blood_flow.unscaled_velocity=(cotd(thetas-90));
mpP.Blood_flow.spread_matrix_adaptive=spread_matrix_adaptive;
mpP.Blood_flow.baseline_var=baseline_var;

objective_microns_per_pixel=user_microns_per_pixel;%
if strcmp(mpP.Header.Scan_Mode,'Line Scan')
    Xfactor=objective_microns_per_pixel/(str2num(mpP.Header.Magnification(1:end-1)));%objective_mag*microns_per_pixel/str2num(mpP.Header.Magnification(2:end));
    Tfactor =mpP.Header.LineRate;
elseif  strcmp(mpP.Header.Scan_Mode,'Arbitrary Scan')
    Xfactor=objective_microns_per_pixel*((mpP.scanData.scanVelocity)/(4/512));%4 volts/512 pixels
    Tfactor=mpP.Header.LineRate;
end

mpP.Xfactor=Xfactor;
mpP.Tfactor=Tfactor;
mpP.Blood_flow.the_decimate=the_decimate;
if isinterleaved==1 %need to divide by 2 because of splitting due to inerleaving
    mpP.Blood_flow.v(1,:)=the_decimate*Tfactor*Xfactor*(cotd(thetas(1,:)-90))/2;
    mpP.Blood_flow.v(2,:)=the_decimate*Tfactor*Xfactor*(cotd(thetas(2,:)-90))/2;
else
    mpP.Blood_flow.v(1,:)=the_decimate*Tfactor*Xfactor*(cotd(thetas(1,:)-90));
    mpP.Blood_flow.v(2,:)=the_decimate*Tfactor*Xfactor*(cotd(thetas(1,:)-90));
end

mpP.Blood_flow.velocity2(1,:)=VelocityCleanUp(mpP.Blood_flow.v(1,:),max_v);

if isinterleaved==1
    mpP.Blood_flow.velocity2(2,:)=VelocityCleanUp(mpP.Blood_flow.v(2,:),max_v);
else
    mpP.Blood_flow.velocity2(2,:)=-(mpP.Blood_flow.velocity2(1,:));
end

v_hold=mpP.Blood_flow.velocity2;
v_out=0.5*(v_hold(1,:)-v_hold(2,:));%take the average of the two directions
mpP.Blood_flow.v_out=v_out;
mpP.Blood_flow.the_t=the_t/(mpP.Header.LineRate);
mpP.Blood_flow.the_Fs=4/window_time_size;
mpP.num_frames

the_date=date;
[fixed_v]=VelocityCleanUpSTDOutliers(mpP.Blood_flow.v_out, 3);% clear outliers > 3standard deviations
mpP.Blood_flow.fixed_v=fixed_v;
params.Fs=1/(mpP.Blood_flow.the_t(2)-mpP.Blood_flow.the_t(1));
params.tapers=[3 5];%to deal with Chronux 2.0
[Shr,thr,fhr]=mtdspecgramc(diff(mpP.Blood_flow.fixed_v),[2 .5],[pi/2],params);

mpP.Blood_flow.velocity_specgram=Shr;
mpP.Blood_flow.velocity_specgram_t=thr;
try %to load raw data
    mpP.Blood_flow.Image=data_cache;
catch
    mpP.Blood_flow.Image=zeros(mpP.num_frames*(mpP.Header.Frame_Height),1);
end



seperability=max(squeeze(mean(mpP.Blood_flow.spread_matrix_adaptive,1))')./mean(squeeze(mean(mpP.Blood_flow.baseline_var,1))');%
mpP.Blood_flow.seperability=seperability; % seperability should be >3
sep_threshold=3;
low_sep=find(seperability<sep_threshold);


save(['mv_mpP_' fname(1:end-4) '_' the_date],'mpP')
end

%% dependent functions
function [fixed_v]=VelocityCleanUpSTDOutliers(v, thresh)
%finds points outside a threshold velocity and interpolated between last
% and next good points, unless at the begining or end, then puts in mean
% velocities
%pd 2-21-08
upper_bound=thresh*std(v)+mean(v);
lower_bound=-thresh*std(v)+mean(v);
badpoints=union(find(v>upper_bound),find(v<lower_bound));
noutliers=length(badpoints);
fixed_v=v;
while noutliers>0
    %badpoints=find(abs(v)>thresh);
    goodpoints=union(setdiff(badpoints+1,badpoints),setdiff(badpoints-1,badpoints));
    for p=1:length(badpoints)
        lastgood=max(find(goodpoints<badpoints(p)));
        nextgood=min(find(goodpoints>badpoints(p)));
        if ((goodpoints(1)>0)&&(goodpoints(end)<=length(fixed_v)))
            if((numel(lastgood)>0)&&(isempty(nextgood)==0))
                fixed_v(badpoints(p))=.5*fixed_v(goodpoints(lastgood))+.5*fixed_v(goodpoints(nextgood));
            elseif((numel(lastgood)>0))
                fixed_v(badpoints(p))=fixed_v(goodpoints(nextgood));
            else
                fixed_v(badpoints(p))=fixed_v(goodpoints(lastgood));
            end
        elseif (goodpoints(1)==0)
            fixed_v(badpoints(p))=mean(fixed_v);%(goodpoints(nextgood));
        elseif (goodpoints(end)==length(fixed_v)+1)
            fixed_v(badpoints(p))=mean(fixed_v);%fixed_v(goodpoints(lastgood));
        elseif numel(goodpoints)==0
            fixed_v(badpoints(p))=-fixed_v;
        end
        
    end
    upper_bound=thresh*std(fixed_v)+mean(fixed_v);
    lower_bound=-thresh*std(fixed_v)+mean(fixed_v);
    badpoints=union(find(fixed_v>upper_bound),find(fixed_v<lower_bound));
    noutliers=length(badpoints);
    
end
end

function [fixed_v]=VelocityCleanUp(v, the_max_v)
%finds points faster than a threshold velocity and interpolated between last
% and next good points
%pd 11-30-07

badpoints=find(abs(v)>the_max_v);
goodpoints=union(setdiff(badpoints+1,badpoints),setdiff(badpoints-1,badpoints));
fixed_v=v;
for p=1:length(badpoints)
    lastgood=max(find(goodpoints<badpoints(p)));
    nextgood=min(find(goodpoints>badpoints(p)));
    if ((goodpoints(1)>0)&&(goodpoints(end)<=length(v)))
        
        if((numel(lastgood)>0)&&(isempty(nextgood)==0))
            fixed_v(badpoints(p))=.5*v(goodpoints(lastgood))+.5*v(goodpoints(nextgood));
        elseif((numel(lastgood)>0))
            fixed_v(badpoints(p))=v(goodpoints(nextgood));
        else
            fixed_v(badpoints(p))=v(goodpoints(lastgood));
        end
    elseif (goodpoints(1)==0)
        fixed_v(badpoints(p))=v(goodpoints(nextgood));
    elseif (goodpoints(end)==length(v)+1)
        fixed_v(badpoints(p))=v(goodpoints(lastgood));
    elseif numel(goodpoints)==0
        fixed_v(badpoints(p))=-v;
    end
    
end
end

function [Header]=ReadMCSTiffInfo(fname)
%extract header information from a MCS file saved as a *tiff using
%mpview tiff exporter
%PJD 5/2017
rawinfo=imfinfo(fname);
the_strings=textscan(rawinfo(1).ImageDescription,'%s','Delimiter', ':');
Header.Filename=rawinfo(1).Filename;
Header.Frame_Count=length(rawinfo);
Header.Scan_Mode=the_strings{1}{17};%
Header.X_Position = the_strings{1}{39};%
Header.Y_Position = the_strings{1}{41};%
Header.Z_Position = the_strings{1}{43};%
Header.Stack_Count = the_strings{1}{29};%
Header.Z_Interval = the_strings{1}{27};%
Header.Averaging_Count = the_strings{1}{31};%
%Header.Repeat_Count = the_strings{1}{20};
Header.Magnification = the_strings{1}{37};%
Header.Rotation = the_strings{1}{35};%
Header.X_Frame_Offset = the_strings{1}{23};%
Header.Y_Frame_Offset = the_strings{1}{25};%
Header.Frame_Duration = the_strings{1}{45};%

Header.Frame_Rate =1/str2num(Header.Frame_Duration(1:end-2));
Header.Filename = rawinfo(1).Filename;
Header.Frame_Width = rawinfo(1).Width;
Header.Frame_Height = rawinfo(1).Height;
Header.Pixel_Clock= round(1/(rawinfo(1).Width*(5/4)*rawinfo(1).Height*(Header.Frame_Rate)*(.05*1e-6))); %not used
Header.Channel_Name1 = 'null';
Header.Channel_Name2 = 'null';
Header.Channel_Name3 = 'null';
Header.Channel_Name4 = 'null';
Header.Enabled1 = 'null';
Header.Enabled2 = 'null';
Header.Enabled3 = 'null';
Header.Enabled4 = 'null';
Header.Input_Range1 = 'null';
Header.Input_Range2 = 'null';
Header.Input_Range3 = 'null';
Header.Input_Range4 = 'null';
Header.Channel_Unit3 = 'null';
Header.Channel_Unit4 = 'null';
Header.Channel_Prefix3 = 'null';
Header.Channel_Prefix4 = 'null';
Header.Conversion_Factor3 = 'null';
Header.Conversion_Factor4 = 'null';
Header.Offset3 = 'null';
Header.Offset4 = 'null';
Header.Data_Point_Per_Frame3 = 'null';
Header.Data_Point_Per_Frame4 = 'null';
end

function [the_image]=LoadTiffConcatenate(the_tiff,the_frames)
% loads the min(frames):maxframes of the_tiff, can conatenates them into
% the_image, a double matrix.  If no frame range given, load all frames.

if isempty(the_frames)
    the_info=imfinfo(the_tiff);
    start_frame=1;%min(the_frames);
    end_frame=length(the_info);
    n_frames=end_frame-start_frame+1;
else
    start_frame=min(the_frames);
    end_frame=max(the_frames);
    n_frames=end_frame-start_frame+1;
end

tiff_file=Tiff(the_tiff);
first_frame=tiff_file.read();
tiff_height=size(first_frame,1);
tiff_width=size(first_frame,2);
the_image=zeros(tiff_height*n_frames,tiff_width);

for n=1:n_frames
    tiff_file.setDirectory(n);
    the_image((1+(n-1)*tiff_height):(n*tiff_height),:)=double(tiff_file.read());
end
tiff_file.close();

end

function [mv_mpP]=linescan_xcov_velocity_TIFF(mv_mpP)
%this function uses the cross correlation method to measure velocity
%using the method of Kim et al, PLoS ONE 2012 DOI:10.1371/journal.pone.0038590
%Patrick Drew 7/2012
%
[~,xc_image,xshift]=linescan_xcov(mv_mpP);
[b,a]=butter(2,200/(2*mv_mpP.Tfactor));%200Hz cutoff
xc_velocity=filtfilt(b,a,xshift*mv_mpP.Xfactor*mv_mpP.Tfactor);% calculate the velocity from the displacement and filter
mv_mpP.Blood_flow.xc_velocity=xc_velocity;%velocity obtained with the cross correlation method
params.Fs=mv_mpP.Tfactor;
params.tapers=[20 39];
[S_xc,f_xc]=mtspectrumc(detrend(xshift),params);
mv_mpP.Blood_flow.S_xcor=S_xc;
mv_mpP.Blood_flow.f_xcor=f_xc;
mv_mpP.Blood_flow.params_xc=params;
mv_mpP.Blood_flow.xc_image=xc_image;

figure(33)
subplot(1,4,1:3)
hold off
imagesc((1:length(mv_mpP.Blood_flow.xc_image))/mv_mpP.Tfactor,(size(mv_mpP.Blood_flow.xc_image,1)/2:-1:-size(mv_mpP.Blood_flow.xc_image,1)/2)*mv_mpP.Tfactor*mv_mpP.Xfactor/1000,mv_mpP.Blood_flow.xc_image)
hold on
plot((1:length(mv_mpP.Blood_flow.xc_velocity))/mv_mpP.Tfactor,mv_mpP.Blood_flow.xc_velocity/1000,'w')
xlabel('time, seconds')
ylabel('velocity, mm/sec')
title('cross-correlation estimate of velocity')
axis xy
subplot(1,4,4)
loglog(f_xc,S_xc)

end

function [maxspot,xc_image,xshift]=linescan_xcov(mv_mpP)
% this calculates the velocity using a cross-correlation technique.  only
% works for very fast moving flow, not capillaries
% maxspot-amplitude of the cross correlation peak.
% xcorr_image- power-spectra normalized cross correlation image
% xshift- shift of the peak away from the center, in pixels
x_spread=round(max(1,.5/mv_mpP.Xfactor));%spatial gaussian with 0.5um std
t_spread=round(mv_mpP.Tfactor/10);%temporal gaussian with 10ms std
the_kernel=gaussian2d(x_spread,t_spread,3*x_spread,3*t_spread);
theimage=double(mv_mpP.Blood_flow.Image');
nlines=size(theimage,2);
npoints=size(theimage,1);
average_line=mean(theimage,2);
xc_image=zeros(size(theimage));
for t=1:nlines
    theimage(:,t)=theimage(:,t)-average_line;%subtract out the background
end
for t=3:nlines
    %take the convolution of the two line-scans normalized by the joint
    %power spectrums
    xc_image(:,t)=ifft(fft(theimage(:,t)).*conj(fft(theimage(:,t-1)))./...
        sqrt(abs(fft(theimage(:,t))).*abs(fft(theimage(:,t-1)))));
end
%convolve with a gaussian space-time kernel to average velocity
xc_image=conv2(fftshift(xc_image,1),the_kernel, 'same');%fftshift puts the lower absolute values for the velocities together in the middel of the matrix,then we filter with a matched kernel
[~,maxspot]=max(xc_image(:,:));%find the peak in the cross correlation
xshift=(round(npoints/2)-maxspot);
end

function [data]=gaussian2d(xstd,ystd,xsize,ysize)
%2-d gaussian kernel for filtering
data=zeros(xsize-1,ysize-1);
x0=xsize/2;
y0=ysize/2;
for x=1:xsize-1
    for y=1:ysize-1
        data(x,y)=exp(-((x-x0)^2)/(2*xstd)-((y-y0)^2)/(2*ystd));
    end
end
data=data/sum(data(:));
end




