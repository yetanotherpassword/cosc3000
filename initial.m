clear all
close all
delimiter = ',';
bomdir="~/bom_bne/";
AllBomData=[];
AllBom=dir (bomdir+"*.csv");
varNames={'Date','MinTemp','MaxTemp','Rainfall','9amTemp','9amRelHumid','9amOktas','9amWindDir','9amWindSpd','9amMSLPres','3pmTemp','3pmRelHumid','3pmOktas','3pmWindDir','3pmWindSpd','3pmMSLPres' };
varTypes = {'char','double','double','double','double',   'double',     'int',     'char',       'int',       'double',    'double',   'double',     'int',     'char',       'int',       'double' } ;
loc=[];
summary_daily_plot=false;
track_plot_daily = false;
summary_monthly_plot=true;
dbgstate=false;;

for j=1:1:size(AllBom,1)
    % details=split(AllBom(j).name,".");
    % dte=string(details(2))
    % yyyy=extractBetween(dte,1,4)
    % mm=extractBetween(dte,5,6)
    opts = detectImportOptions(bomdir+AllBom(j).name,'Delimiter',delimiter,'PartialFieldRule','fill');
    opts.SelectedVariableNames = opts.SelectedVariableNames([2, 3, 4, 5, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]);
    T = readtable(bomdir+AllBom(j).name,opts);
    T.Properties.VariableNames=varNames;
    AllBomData=[AllBomData; T];
    fid = fopen(bomdir+AllBom(j).name,'rt');
    summary1 = fgetl(fid);
    fclose(fid);
    loc=[loc; string(summary1)];
end



AllDays=dir ("~/bn_data/");



altcolr=[0 1 0];
varNames = {'GeoHeight', 'Elevation', 'RFS_Detection_Time','RFS_Range','RFS_Azimuth', 'BaroAlt', 'RFS_Velocity','RFS_Heading','RangeError','App'} ;
varTypes = {'double','double','double','double','double','double','double','double','double','char'} ;

dataStartLine = 1;
extraColRule = 'ignore';
minalt=16000;
maxalt=0;
colbar=[];
altbar=[];

sampletracks_3pm=0;
sampletracks_9am=0;

extrap_rng_err_9am=[];
extrap_rng_err_3pm=[];

y_3pm=[];
ym_3pm=[];
ys_3pm=[];
y_9am=[];
ym_9am=[];
ys_9am=[];
actual_T1_samples=[];
actual_T1_Range=[];
actual_T1_RangeErr=[];
ytraps_3pm=[];
ytraps_9am=[];
cols_3pm=[];
alts_3pm=[];
cols_9am=[];
alts_9am=[];

colour_approach_9am=[];
actual_T1_samples_9am=[];
actual_T1_Range_9am=[];
actual_T1_RangeErr_9am=[];
colour_approach_3pm=[];
actual_T1_samples_3pm=[];
actual_T1_Range_3pm=[];
actual_T1_RangeErr_3pm=[];


daily_9am=[];
daily_9am_day=[];
daily_9am_temp=[];
daily_9am_rh=[];
daily_3pm=[];
daily_3pm_day=[];
daily_3pm_temp=[];
daily_3pm_rh=[];

monthly_3pm=[];
monthly_3pm_day=[];
monthly_3pm_temp=[];
monthly_3pm_rh=[];

monthly_9am=[];
monthly_9am_day=[];
monthly_9am_temp=[];
monthly_9am_rh=[];

yearly_3pm=[];
yearly_3pm_day=[];
yearly_3pm_temp=[];
yearly_3pm_rh=[];

yearly_9am=[];
yearly_9am_day=[];
yearly_9am_temp=[];
yearly_9am_rh=[];

ac_stats=["ModeSId","Reports", "DeltaErr", "MinAlt", "MaxAlt", "MinRng", "MaxRng", "MinElev", "MaxElev", "Start", "Stop"];
for i = 1:1:size(AllDays,1)
    if length(AllDays(i).name) ~= 10
        continue
    else
        if (AllDays(i).isdir)
            ThisDate=sprintf ("%s",AllDays(i).name);
            fprintf("%%%%%%%%%% Doing %s\n",ThisDate);
            CurrPath="~/bn_data/"+ThisDate+"/SIC_045_BNE/TGT/CSV/";
            AllTracks=dir (CurrPath+"SIC*.csv");
            NumberOfTracksThisDay=size(AllTracks,1);
            sampletracks_9am=0;
            sampletracks_3pm=0;
            for j=1:1:NumberOfTracksThisDay
                if length(AllTracks(j).name)>2
                    % (>2 means: ignore current and parent dir filenames)
                    ThisTrackFname=CurrPath+AllTracks(j).name;
                    details=split(AllTracks(j).name,"_");
                    qos=is_qos(details);

                    msid=string(details(6));
                    opts = detectImportOptions(ThisTrackFname,'Delimiter',delimiter,'PartialFieldRule','keep');
                    opts.SelectedVariableNames = opts.SelectedVariableNames([9, 27, 29, 31, 32, 34, 38, 39, 44, 79]);

                    % Read the stat summary -extra fields on 1st data line
                    fid = fopen(ThisTrackFname,'rt');
                    summary1 = fscanf(fid, '%s\n');
                    fclose(fid);
                    a=split(summary1,',');
                    RangeBias=str2double(a(163)); % ie field 71 on 2nd line


                    T = readtable(ThisTrackFname,opts);
                    T.Properties.VariableNames=varNames;
                    today=posixtime(datetime(ThisDate));
                    todaylessday=today-86400;
                    todayplusday=today+86400;
                    % Do yesterday first as will then continue to today
                    % incase it is available
                    yester_dat=datetime(todaylessday, 'convertfrom', 'posixtime', 'Format', 'yyyy-MM-dd');
                    yester_str=string(yester_dat);
                    tommor_dat=datetime(todayplusday, 'convertfrom', 'posixtime', 'Format', 'yyyy-MM-dd');
                    tomor_str=string(tommor_dat);
                    tomor_mth=extractBetween(tomor_str,1,7);
                    today_mth=extractBetween(ThisDate,1,7);
                    yesterday_3pm=getTime(5*3600, T);
                    today_9am=getTime(23*3600, T);
                    for q=1:1:2
                        if (q == 1)
                            this_sample=yester_str+" 3pm "+msid;
                            if ~isempty(yesterday_3pm)
                                T1=T(yesterday_3pm,:);
                            else
                                debug_out(dbgstate, sprintf("Skipping %s as no samples in this time\n",this_sample));
                                continue
                            end
                        else
                            this_sample=ThisDate+" 9am "+msid;
                            if ~isempty(today_9am)
                                T1=T(today_9am,:);
                            else
                                debug_out(dbgstate, sprintf("Skipping %s as no samples in this time\n",this_sample));
                                continue
                            end
                        end


                        samples=size(T1.RFS_Range,1);
                        [max_range, max_idx]=max(T1.RFS_Range);
                        [min_range, min_idx]=min(T1.RFS_Range);
                        delta_error=abs(diff(T1.RangeError));
                        max_delta_error = max(delta_error);
                        if (samples < 200 || max_range-min_range < 200 || max_delta_error > 400)
                            debug_out(dbgstate, sprintf("Rejecting %s as samples=%d, delta range=%f max delta error=%f\n",this_sample,samples,max_range-min_range,max_delta_error));
                            continue;
                        end

                        yesbomref=find(AllBomData.Date==yester_dat);
                        if isempty(yesbomref)
                            debug_out(dbgstate, sprintf("Error missing %s in BOM Data skipping...\n",this_sample));
                            continue;
                        end
                        bomref=find(AllBomData.Date==ThisDate);
                        if isempty(bomref)
                            debug_out(dbgstate, sprintf("Error missing  %s in BOM Data\n",this_sample));
                            continue;
                        end


                        debug_out(dbgstate, sprintf("\n******** PASSED checks: Proceeding with %s (Number %d of %d)\n*********\n",this_sample,j,NumberOfTracksThisDay));


                        minalt=min(minalt,min(T1.GeoHeight));
                        maxalt=max(maxalt,max(T1.GeoHeight));
                        mintime=min(T1.RFS_Detection_Time);
                        maxtime=max(T1.RFS_Detection_Time);
                        minelev=min(T1.Elevation);
                        maxelev=max(T1.Elevation);

                        ac_data=[msid samples max_delta_error minalt maxalt  min_range max_range minelev maxelev mintime maxtime];

                        %Get Elevation Data (Most interested in <5 and <2 deg
                        %col1=[0.9290 0.6940 0.1250]; % Gold >10 deg elev
                        col1=[1 228/255 149/255];
                        col2=[1 0.5 0.5]; % Red between 5 and 10 deg elev
                        %col3=[0.4660 0.6740 0.1880]; % Green between 2 and 5 deg elevation
                        col3=[7/255 204/255 7/255];
                        col4=[0 0 1]; % Blue less than 2 deg elevation
                        col_elev=zeros(samples,3);
                        col_elev(:,:)=col_elev(:,:)+col1;
                        col_lt10=find(T1.Elevation<=10);
                        col_lt5=find(T1.Elevation<=5);
                        col_lt2=find(T1.Elevation<=2);
                        col_elev(col_lt10,:)=zeros(size(col_lt10,1),3)+col2;
                        col_elev(col_lt5,:)=zeros(size(col_lt5,1),3)+col3; % 1==OUT Bound, 0==IN bound
                        col_elev(col_lt2,:)=zeros(size(col_lt2,1),3)+col4;



                        % TODO: Reflection : Pain for 3pm being 'yesterday' as UTC time used
                        %                    Problem with sample size varying, cant put into matlab
                        %                    matrix unless all same size (padding) so encoding with
                        %                    in col vector using a size vertor

                        if track_plot_daily == true
                            f=figure;
                            f.Renderer = 'painters';
                            p1=plot(T1.RFS_Range,T1.RangeError-RangeBias,':k','LineWidth', 1);
                            hold on
                        end
                        coefficients = polyfit(T1.RFS_Range,T1.RangeError-RangeBias, 2);
                        xFit = linspace(min_range, max_range, samples);
                        yFit = polyval(coefficients , xFit);
                        txt1=sprintf("%s\n",ac_stats);
                        txt2=sprintf("%s\n",ac_data);
                        fv = [1:samples-1;2:samples]';
                        if track_plot_daily == true
                            p2=patch('faces',fv,'vertices',[xFit; yFit]',...
                                'faceVertexCData',col_elev,...
                                'edgecolor','flat',...
                                'linewidth',2)

                            %scatter(xFit, yFit,1,col_elev, 'LineWidth', 2);
                            hold on;


                            colormap([col1; col2; col3; col4]);


                            hcb=colorbar('Ticks',[0.125,0.375,0.625,0.875],...
                                'TickLabels',{'>10^o','<=10^o','<=5^o','<=2^o'})

                            colorTitleHandle = get(hcb,'Title');
                            titleString = sprintf('Aircraft Elevation Angle\n(to Radar)');
                            set(colorTitleHandle ,'String',titleString);
                        end
                        titlename=sprintf("Range vs Range Error (Actual and Regession)\n");
                        if q==1
                            Temp_3pm=AllBomData{yesbomref,11};
                            RelHum_3pm=AllBomData{yesbomref,12};
                            Date_3pm=yester_str;
                            Period_3pm = sprintf("3pm +/- 1 hr %s      %02.1f^oC      %d%% Rel Humid",Date_3pm,Temp_3pm,RelHum_3pm);

                            %Period="3pm "+string(yester)+" "+string(AllBomData{yesbomref,11})+"^oC "+string(AllBomData{yesbomref,12})+"% Rel Humd";
                            actual_T1_samples_3pm=[actual_T1_samples_3pm; samples];
                            actual_T1_Range_3pm=[actual_T1_Range_3pm; T1.RFS_Range];
                            actual_T1_RangeErr_3pm=[actual_T1_RangeErr_3pm; T1.RangeError-RangeBias];
                            fname="~/"+msid+"_3pm_"+yester_str+".png";

                        else
                            Temp_9am=AllBomData{bomref,5};
                            RelHum_9am=AllBomData{bomref,6};
                            Date_9am=ThisDate;
                            Period_9am = sprintf("9am +/- 1 hr %s      %02.1f^oC      %d%% Rel Humid",Date_9am,Temp_9am,RelHum_9am);


                            %Period = "9am "+ThisDate+" "+string(AllBomData{bomref,5})+"^oC "+string(AllBomData{bomref,6})+"% Rel Humd";
                            actual_T1_samples_9am=[actual_T1_samples_9am; samples];
                            actual_T1_Range_9am=[actual_T1_Range_9am; T1.RFS_Range];
                            actual_T1_RangeErr_9am=[actual_T1_RangeErr_9am; T1.RangeError-RangeBias];
                            fname="~/"+msid+"_9am_"+ThisDate+".png";
                        end




                        if track_plot_daily == true
                            xlabel("Range to Aircraft (Nautical Miles)");
                            ylabel("Range Error (metres)");
                            legend("Actual Error","Regression Error");
                            if q==1
                                title(titlename+Period_3pm)
                            else
                                title(titlename+Period_9am)
                            end
                            yl=ylim;
                            text(220,yl(2)-0.2*(yl(2)-yl(1)),txt1)
                            text(235,yl(2)-0.2*(yl(2)-yl(1)),txt2)
                            box on
                            set(gcf, 'Position', get(0, 'Screensize'));
                            saveas(gcf,fname);
                            close
                        end

                        % Now extrapolate the range error over the whole radar range
                        xtrap=[0:1:250];

                        ytrap = polyval(coefficients, xtrap);
                        % Save these interpolated errors at sample points for later comparison




                        %extrap_rng_err(find(isnan(extrap_rng_err)))=[];


                        [rerr_avg rerr_std rerr_median] = stats(ytrap);


                        if q==1
                            y_3pm=[y_3pm; rerr_avg ];
                            ym_3pm=[ym_3pm; rerr_median];
                            ys_3pm=[ys_3pm; rerr_std];
                            extrap_rng_err_3pm=[extrap_rng_err_3pm; ytrap];
                            sampletracks_3pm = sampletracks_3pm + 1;
                        else
                            y_9am=[y_9am; rerr_avg];
                            ym_9am=[ym_9am; rerr_median ];
                            ys_9am=[ys_9am; rerr_std];
                            extrap_rng_err_9am=[extrap_rng_err_9am; ytrap];
                            sampletracks_9am = sampletracks_9am + 1;
                        end
                    end

                end
            end
            %        figure;
            %
            %                    coefficients = polyfit(actual_T1_Range_9am,actual_T1_RangeErr_9am, 2);
            %                    xFit = linspace(0, 250, size(actual_T1_Range_9am,1));
            %                    yFit = polyval(coefficients , xFit);
            %                    txt=sprintf("9am : %d Tracks for %s",sampletracks_9am,ThisDate);
            %                    scatter(xFit, yFit,20);
            %                    hold on;
            %                    title("Range VS Interpolated Range Error (Incoming/Outgoing)"+txt)
            %                    xlabel("Range (Nautical Miles)");
            %                    ylabel("Range Error (metres) (Interpolated)")
            %                    inoutcolr=[0 1 0; 0 1 1]; colormap(inoutcolr); colorbar
            %                    grid on
            if summary_daily_plot
                plot_stats(extrap_rng_err_3pm,"Daily 3pm"+yester_str);
                plot_stats(extrap_rng_err_9am,"Daily 3pm"+today_str);
                %[y2 y2s y2m] = stats(extrap_rng_err_3pm); % check if same or not
            end
            daily_3pm=[daily_3pm; extrap_rng_err_3pm];
            daily_3pm_day=[daily_3pm_day; Date_3pm];
            daily_3pm_temp=[daily_3pm_temp; Temp_3pm]
            daily_3pm_rh=[daily_3pm_rh; RelHum_3pm]

            daily_9am=[daily_9am; extrap_rng_err_9am];
            daily_9am_day=[daily_9am_day; Date_9am];
            daily_9am_temp=[daily_9am_temp; Temp_9am]
            daily_9am_rh=[daily_9am_rh; RelHum_9am]



        else
            debug_out(dbgstate, sprintf("%s is not a dir, skipping...\n",AllDays(i).name));
            continue
        end
    end
    % Get monthly max and min temp and humidity
    [max9t maxidx9t]=max(daily_9am_temp);
    [min9t minidx9t]=min(daily_9am_temp);
    [max3t maxidx3t]=max(daily_3pm_temp);
    [min3t minidx3t]=min(daily_3pm_temp);

    [max9h maxidx9h]=max(daily_9am_rh);
    [min9h minidx9h]=min(daily_9am_rh);
    [max3h maxidx3h]=max(daily_3pm_rh);
    [min3h minidx3h]=min(daily_3pm_rh);

    % plot these
    if summary_monthly_plot
        if tomor_mth ~= today_mth
            plot_stats(daily_9am,"Monthly 9am "+today_mth);
            plot_stats(daily_3pm,"Monthly 3pm"+today_mth); % use 9am as we will use that month (next 3pm will be used for next month)

            yearly_3pm=[yearly_3pm; monthly_3pm];
            yearly_3pm_day=[yearly_3pm_day; monthly_3pm_day];
            yearly_3pm_temp=[yearly_3pm_temp; monthly_3pm_temp];
            yearly_3pm_rh=[yearly3_pm_rh; monthly_3pm_rh];

            yearly_9am=[yearly_9am; monthly_9am];
            yearly_9am_day=[yearly_9am_day; monthly_9am_day];
            yearly_9am_temp=[yearly_9am_temp; monthly_9am_temp];
            yearly_9am_rh=[yearly_9am_rh; monthly_9am_rh];

            monthly_3pm=[];
            monthly_3pm_day=[];
            monthly_3pm_temp=[];
            monthly_3pm_rh=[];

            monthly_9am=[];
            monthly_9am_day=[];
            monthly_9am_temp=[];
            monthly_9am_rh=[];

        end

    end
    % save yearly for overall average to compare later
    monthly_3pm=[monthly_3pm; daily_3pm];
    monthly_3pm_day=[monthly_3pm_day; daily_3pm_day];
    monthly_3pm_temp=[monthly_3pm_temp; daily_3pm_temp];
    monthly_3pm_rh=[monthly_3pm_rh; daily_3pm_rh];

    monthly_9am=[monthly_9am; daily_9am];
    monthly_9am_day=[monthly_9am_day; daily_9am_day];
    monthly_9am_temp=[monthly_9am_temp; daily_9am_temp];
    monthly_9am_rh=[monthly_9am_rh; daily_9am_rh];

    daily_9am=[];
    daily_9am_day=[];
    daily_9am_temp=[];
    daily_9am_rh=[];

    daily_3pm=[];
    daily_3pm_day=[];
    daily_3pm_temp=[];
    daily_3pm_rh=[];
end
%plot(x,y, '-r');
%plot(x,ym, '-b');

%legend("3pm Average of population at Range X","Median of population at Range X");

function [avg stdev median1]=stats(array)
avg=mean(array);
stdev=std(array);
median1=median(array);
end

function  colvect=colorAtAlt(i)
newalt = i*65535/15000;
c1=newalt/256;
c2 = rem(newalt,256);
colvect=[1 c2/256 c1/256];
end

function [col avgalt]=getAlt(i, T)
[TR501 fiftyidx]=find(T.RFS_Range>i-0.1);
if isempty(fiftyidx)
    col=[1 1 1];
    avgalt=-1;
else
    [TR502 fiftyidx2]=find(TR501 < i+0.1);
    if isempty(fiftyidx2)
        col = [1 1 1];
        avgalt=-1;
    else
        alts1 = T.GeoHeight(fiftyidx);
        alts2 = alts1(fiftyidx2);
        avgalt=mean(alts2);
        col=colorAtAlt(avgalt);
    end
end
end

function timeidx=getTime(i, T)
if (i-3600 >= 0 && i+3600 <= 86400)
    j1=i-3600;
    j2=i+3600;
    timeidx0=find(T.RFS_Detection_Time>j1);
    if isempty(timeidx0)
        timeidx=find(T.RFS_Detection_Time<j2);
    else
        tmp=ones(size(T,1),1)*90000;
        tmp(timeidx0)=T.RFS_Detection_Time(timeidx0);
        timeidx=find(tmp<j2);
    end
end
end

function debug_out(doit, str)
if doit
    fprintf("****** "+str);
end
end

function plot_stats(errs, titnam)
figure
len=size(errs,2)
x=[0:1:250];
y=mean(errs,'omitnan');
ym=median(errs,'omitnan');
ys=std(errs,'omitnan');
ymax=max(errs,[],'omitnan');
ymin=min(errs,[],'omitnan');
errorbar(x,y,ys);
hold on
errorbar(x,ym,ymin,ymax)
legend("Average (+/- STD)","Median (Max,Min)")
title(titnam) % first to last of daily_9am_day
%errorbar(x,ym,[rerr_0_std rerr_50_std rerr_100_std rerr_150_std rerr_200_std rerr_250_std ]);
%[y2 y2s y2m] = stats(daily9am); % check if same or not
end

function q=is_qos(d)
if length(d)>8
    q="yes";
else
    q="no";
end
end