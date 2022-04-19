clear all
close all
delimiter = ',';
bomdir="~/bom_bne/";
AllBomData=[];
AllBom=dir (bomdir+"*.csv");
varNames={'Date','MinTemp','MaxTemp','Rainfall','9amTemp','9amRelHumid','9amOktas','9amWindDir','9amWindSpd','9amMSLPres','3pmTemp','3pmRelHumid','3pmOktas','3pmWindDir','3pmWindSpd','3pmMSLPres' };
varTypes = {'char','double','double','double','double',   'double',     'int',     'char',       'int',       'double',    'double',   'double',     'int',     'char',       'int',       'double' } ;
loc=[];

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



AllDays=dir ("~/bn_data/")



altcolr=[0 1 0]
varNames = {'GeoHeight', 'Elevation', 'RFS_Detection_Time','RFS_Range','RFS_Azimuth', 'BaroAlt', 'RFS_Velocity','RFS_Heading','RangeError','App'} ;
varTypes = {'double','double','double','double','double','double','double','double','double','char'} ;

dataStartLine = 1;
extraColRule = 'ignore';
minalt=16000;
maxalt=0;
colbar=[];
altbar=[];
sampletracks=0;
sampletracks_3pm=0;
sampletracks_9am=0;

extrap_rng_err_0=[];
extrap_rng_err_50=[];
extrap_rng_err_100=[];
extrap_rng_err_150=[];
extrap_rng_err_200=[];
extrap_rng_err_250=[];
       y_3pm=[];
ym_3pm=[];
y_9am=[];
ym_9am=[];
actual_T1_samples=[];
actual_T1_Range=[];
actual_T1_RangeErr=[];
ytraps_3pm=[];
ytraps_9am=[];
cols_3pm=[];
alts_3pm=[];
cols_9am=[];
alts_9am=[];
Period_9am=[];
Period_3pm=[];
colour_approach_9am=[];  
actual_T1_samples_9am=[];
actual_T1_Range_9am=[];
actual_T1_RangeErr_9am=[];
colour_approach_3pm=[];  
actual_T1_samples_3pm=[];
actual_T1_Range_3pm=[];
actual_T1_RangeErr_3pm=[];
for i = 1:1:size(AllDays,1)

    if (AllDays(i).isdir && length(AllDays(i).name) == 10 && AllDays(i).name=="2021-06-06")
       ThisDate=sprintf ("%s",AllDays(i).name);

       CurrPath="~/bn_data/"+ThisDate+"/SIC_045_BNE/TGT/CSV/";
       AllTracks=dir (CurrPath+"SIC*.csv");
       NumberOfTracksThisDay=size(AllTracks,1);

       for j=1:1:NumberOfTracksThisDay
           if length(AllTracks(j).name)>2
               % (>2 means: ignore current and parent dir filenames)
               ThisTrackFname=CurrPath+AllTracks(j).name;
               details=split(AllTracks(j).name,"_");
               if length(details)>8
                   qos="yes";
               else
                   qos="no";
               end
               msid=string(details(6));
               opts = detectImportOptions(ThisTrackFname,'Delimiter',delimiter,'PartialFieldRule','keep');
               opts.SelectedVariableNames = opts.SelectedVariableNames([9, 27, 29, 31, 32, 34, 38, 39, 44, 79]);  

               % Read the stat summary (extra fields only on first data
               % line)
               fid = fopen(ThisTrackFname,'rt');
               summary1 = fscanf(fid, '%s\n');
               fclose(fid);
               a=split(summary1,',');
               RangeBias=str2double(a(163)); % ie field 71 on 2nd line

               
               T = readtable(ThisTrackFname,opts);
               T.Properties.VariableNames=varNames;

               yesterday_3pm=getTime(5*3600, T);
               today_9am=getTime(23*3600, T);
               for q=1:1:2
                   if (q == 1)
                       if ~isempty(yesterday_3pm)
                           T1=T(yesterday_3pm,:);
                       else
                           continue
                       end
                   else
                       if ~isempty(today_9am)
                           T1=T(today_9am,:);
                       else
                           continue
                       end
                   end

                   samples=size(T1.RFS_Range,1);
                   [max_range, max_idx]=max(T1.RFS_Range);
                   [min_range, min_idx]=min(T1.RFS_Range);
                   delta_error=abs(diff(T1.RangeError));
                   max_delta_error = max(delta_error);
                   if (samples < 200 || max_range-min_range < 200 || max_delta_error > 400)
                       continue;
                   end
                   today=posixtime(datetime(ThisDate));
                   todaylessday=today-86400;
                   % Do yesterday first as will then continue to today
                   % incase it is available
                   yester=datetime(todaylessday, 'convertfrom', 'posixtime', 'Format', 'yyyy-MM-dd');
                   yesbomref=find(AllBomData.Date==yester);
                   if isempty(yesbomref)
                        fprintf("Error missing yesterday %s in BOM Data skipping...\n",yester);
                        continue;
                   end
                   bomref=find(AllBomData.Date==ThisDate);
                   if isempty(bomref)
                       fprintf("Error missing day %s in BOM Data\n",ThisDate)
                       continue;
                   end

                   minalt=min(minalt,min(T1.GeoHeight));
                   maxalt=max(maxalt,max(T1.GeoHeight));
                   mintime=min(T1.RFS_Detection_Time);
                   maxtime=max(T1.RFS_Detection_Time);
                   minelev=min(T1.Elevation);
                   maxelev=max(T1.Elevation);
                   ac_stats=["ModeSId","Reports", "DeltaErr", "MinAlt", "MaxAlt", "MinRng", "MaxRng", "MinElev", "MaxElev", "Start", "Stop"]
                   ac_data=[msid samples max_delta_error minalt maxalt  min_range max_range minelev maxelev mintime maxtime];
%col1=[0.9290 0.6940 0.1250]; % Gold >10 deg elev
col1=[1 228/255 149/255]
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

 %                  colour_approach=[zeros(samples,1) ones(samples,1) col_approach];

%                    figure(1);
%                    plot(T1.RFS_Range,T1.RangeError,'LineWidth', 2);
%                    %scatter(T.RFS_Range,T.RangeError,1,colr, 'LineWidth', 2);
%                    hold on;

% TODO: Reflection : Pain for 3pm being 'yesterday' as UTC time used
%                    Problem with sample size varying, cant put into matlab
%                    matrix unless all same size (padding) so encoding with
%                    in col vector using a size vertor

                    f=figure;
                    f.Renderer = 'painters';
                    p1=plot(T1.RFS_Range,T1.RangeError-RangeBias,':k','LineWidth', 1);
                    hold on
                    coefficients = polyfit(T1.RFS_Range,T1.RangeError-RangeBias, 2);
                    xFit = linspace(min_range, max_range, samples);
                    yFit = polyval(coefficients , xFit);
                    txt1=sprintf("%s\n",ac_stats);
                    txt2=sprintf("%s\n",ac_data);
                    fv = [1:samples-1;2:samples]';
p2=patch('faces',fv,'vertices',[xFit; yFit]',...
'faceVertexCData',col_elev,...
'edgecolor','flat',...
'linewidth',2)

                    %scatter(xFit, yFit,1,col_elev, 'LineWidth', 2);
                    hold on;


                    %colormap([0.9290 0.6940 0.1250; 1 0 0; 0.4660 0.6740 0.1880; 0 0 1]);
                    colormap([col1; col2; col3; col4]);


                   hcb=colorbar('Ticks',[0.125,0.375,0.625,0.875],...
                            'TickLabels',{'>10^o','<=10^o','<=5^o','<=2^o'})

             colorTitleHandle = get(hcb,'Title');
titleString = sprintf('Aircraft Elevation Angle\n(to Radar)');
set(colorTitleHandle ,'String',titleString);

                    if q==1
                       Period = sprintf("3pm +/- 1 hr %s      %s^oC      %s%% Rel Humid",string(yester),string(AllBomData{yesbomref,11}),string(AllBomData{yesbomref,12}));                    

                       %Period="3pm "+string(yester)+" "+string(AllBomData{yesbomref,11})+"^oC "+string(AllBomData{yesbomref,12})+"% Rel Humd";
                      actual_T1_samples_3pm=[actual_T1_samples_3pm; samples];
                       actual_T1_Range_3pm=[actual_T1_Range_3pm; T1.RFS_Range];
                       actual_T1_RangeErr_3pm=[actual_T1_RangeErr_3pm; T1.RangeError-RangeBias];
                       fname="~/"+msid+"_3pm_"+string(yester)+".png";

                    else  
                        Period = sprintf("9am +/- 1 hr %s      %s^oC      %s%% Rel Humid",ThisDate,string(AllBomData{bomref,5}),string(AllBomData{bomref,6}));                    

                       %Period = "9am "+ThisDate+" "+string(AllBomData{bomref,5})+"^oC "+string(AllBomData{bomref,6})+"% Rel Humd";                        
                        actual_T1_samples_9am=[actual_T1_samples_9am; samples];
                       actual_T1_Range_9am=[actual_T1_Range_9am; T1.RFS_Range];
                       actual_T1_RangeErr_9am=[actual_T1_RangeErr_9am; T1.RangeError-RangeBias];
                       fname="~/"+msid+"_9am_"+ThisDate+".png";
                    end

                                   

                    titlename=sprintf("Range vs Range Error (Actual and Regession)\n");
                    title(titlename+Period);
                    xlabel("Range to Aircraft (Nautical Miles)");
                    ylabel("Range Error (metres)");
                    legend("Actual Error","Regression Error")
                                        yl=ylim;
                    text(220,yl(2)-0.2*(yl(2)-yl(1)),txt1)
                    text(235,yl(2)-0.2*(yl(2)-yl(1)),txt2)
                    box on
                    set(gcf, 'Position', get(0, 'Screensize'));
                    saveas(gcf,fname);
                    close


%                    figure(3);

                    [colA altA]=getAlt(5, T1);
                    [colB altB]=getAlt(50, T1);
                    [colC altC]=getAlt(100, T1);
                    [colD altD]=getAlt(150, T1);
                    [colE altE]=getAlt(200, T1);
                    [colF altF]=getAlt(250, T1);
                    colbar=[colbar; colA; colB; colC; colD; colE; colF];
                    altbar=[altbar; altA; altB; altC; altD; altE; altF];

% Take sample ranges for error comparison
                   xtrap = [0 50 100 150 200 250];
% Interpolate quadratic line of best fit
                   ytrap = polyval(coefficients, xtrap);
% Save these interpolated errors at sample points for later comparison
                   extrap_rng_err_0=[extrap_rng_err_0; ytrap(1)];
                   extrap_rng_err_50=[extrap_rng_err_50; ytrap(2)];
                   extrap_rng_err_100=[extrap_rng_err_100; ytrap(3)];
                   extrap_rng_err_150=[extrap_rng_err_150; ytrap(4)];
                   extrap_rng_err_200=[extrap_rng_err_200; ytrap(5)];
                   extrap_rng_err_250=[extrap_rng_err_250; ytrap(6)];

                   %scatter(xtrap, ytrap, 20, [colA; colB; colC; colD; colE; colF],"filled");
                   %hold on;



                   if (q == 1)
                        Period="3pm "+string(yester)+" "+string(AllBomData{yesbomref,11})+"c "+string(AllBomData{yesbomref,12})+"% ";
                        Period_3pm=[Period_3pm; Period];
                        ytraps_3pm=[ytraps_3pm; ytrap Period];
                        cols_3pm=[cols_3pm; colA; colB; colC; colD; colE; colF];
                        alts_3pm=[alts_3pm; altA; altB; altC; altD; altE; altF];
                        sampletracks_3pm=sampletracks_3pm+1;
                   else
                       Period = "9am "+ThisDate+" "+string(AllBomData{bomref,5})+"c "+string(AllBomData{bomref,6})+"% ";
                       Period_9am=[Period_9am; Period];
                       ytraps_9am=[ytraps_9am; ytrap ];
                       cols_9am=[cols_9am; colA; colB; colC; colD; colE; colF];
                       alts_9am=[alts_9am; altA; altB; altC; altD; altE; altF];
                       sampletracks_9am=sampletracks_9am+1;
                   end
                   ttrk=sprintf("%s (+/- 1 Hr)\nSample %d Aircraft\n",Period,sampletracks);
%                   figure(1)
%                    grid on;
%                    title("Range VS Range Error : "+ttrk)
%                    xlabel("Range (Nautical Miles)");
%                    ylabel("Range Error (metres)")
%                    %inoutcolr=[0 1 0; 0 1 1];
%                    %colormap(inoutcolr);
%                    %colorbar

%                    figure(2)
%                    title("Range VS Interpolated Range Error (Incoming/Outgoing)"+ttrk)
%                    xlabel("Range (Nautical Miles)");
%                    ylabel("Range Error (metres) (Interpolated)")
%                    inoutcolr=[0 1 0; 0 1 1];
%                    colormap(inoutcolr);
%                    colorbar
%                    grid on

%                    figure(3);
%                    title("Population of Range VS Interpolated Range Error at Altitude"+ttrk)
%                    xlabel("Range (Nautical Miles)");
%                    ylabel("Range Error (metres) (Interpolated)")
%                    %altcolor=[[1[minalt:1:maxalt][minalt:1:maxalt]];
%                    [saltbar i]=sort(altbar)
%                    newcolbar=colbar(i,:)
%                    colormap(newcolbar);
%                    colorbar
%                    set(gca,'clim',saltbar([1,end]))

                   extrap_rng_err_0(find(isnan(extrap_rng_err_0)))=[];
                   extrap_rng_err_50(find(isnan(extrap_rng_err_50)))=[];
                   extrap_rng_err_100(find(isnan(extrap_rng_err_100)))=[];
                   extrap_rng_err_150(find(isnan(extrap_rng_err_150)))=[];
                   extrap_rng_err_200(find(isnan(extrap_rng_err_200)))=[];
                   extrap_rng_err_250(find(isnan(extrap_rng_err_250)))=[];

                   [rerr_0_avg rerr_0_std rerr_0_median] = stats(extrap_rng_err_0);
                   [rerr_50_avg rerr_50_std rerr_50_median] = stats(extrap_rng_err_50);
                   [rerr_100_avg rerr_100_std rerr_100_median] = stats(extrap_rng_err_100);
                   [rerr_150_avg rerr_150_std rerr_150_median] = stats(extrap_rng_err_150);
                   [rerr_200_avg rerr_200_std rerr_200_median] = stats(extrap_rng_err_200);
                   [rerr_250_avg rerr_250_std rerr_250_median] = stats(extrap_rng_err_250);
                   if q==1
                       y_3pm=[y_3pm; rerr_0_avg rerr_50_avg rerr_100_avg rerr_150_avg rerr_200_avg rerr_250_avg];
                       ym_3pm=[ym_3pm; rerr_0_median rerr_50_median rerr_100_median rerr_150_median rerr_200_median rerr_250_median];
                   else
                       y_9am=[y_9am; rerr_0_avg rerr_50_avg rerr_100_avg rerr_150_avg rerr_200_avg rerr_250_avg];
                       ym_9am=[ym_9am; rerr_0_median rerr_50_median rerr_100_median rerr_150_median rerr_200_median rerr_250_median];
                   end         
               end

           end
       end
       figure;

                   coefficients = polyfit(actual_T1_Range_9am,actual_T1_RangeErr_9am, 2);
                   xFit = linspace(0, 250, size(actual_T1_Range_9am,1));
                   yFit = polyval(coefficients , xFit);
                   txt=sprintf("9am : %d Tracks for %s",sampletracks_9am,ThisDate);
                   scatter(xFit, yFit,20);
                   hold on;
                   title("Range VS Interpolated Range Error (Incoming/Outgoing)"+txt)
                   xlabel("Range (Nautical Miles)");
                   ylabel("Range Error (metres) (Interpolated)")
                   %inoutcolr=[0 1 0; 0 1 1];
                   %colormap(inoutcolr);
                   %colorbar
                   grid on
    end
  

end
                figure
                   x=[0 50 100 150 200 250];
                   y=[rerr_0_avg rerr_50_avg rerr_100_avg rerr_150_avg rerr_200_avg rerr_250_avg];
                   ym=[rerr_0_median rerr_50_median rerr_100_median rerr_150_median rerr_200_median rerr_250_median];
                   errorbar(x,y,[rerr_0_std rerr_50_std rerr_100_std rerr_150_std rerr_200_std rerr_250_std ]);
                   hold on
                   %errorbar(x,ym,[rerr_0_std rerr_50_std rerr_100_std rerr_150_std rerr_200_std rerr_250_std ]);

                   %plot(x,y, '-r');
                   plot(x,ym, '-b');
                   legend("Average of population at Range X","Median of population at Range X");

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

