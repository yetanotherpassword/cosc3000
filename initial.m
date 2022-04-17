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
varNames = {'GeoHeight', 'RFS_Detection_Time','RFS_Range','RFS_Azimuth', 'BaroAlt', 'RFS_Velocity','RFS_Heading','RangeError','App'} ;
varTypes = {'double','double','double','double','double','double','double','double','char'} ;

dataStartLine = 1;
extraColRule = 'ignore';
minalt=16000;
maxalt=0;
colbar=[];
altbar=[];
sampletracks=0;

extrap_rng_err_0=[];
extrap_rng_err_50=[];
extrap_rng_err_100=[];
extrap_rng_err_150=[];
extrap_rng_err_200=[];
extrap_rng_err_250=[];

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
               opts.SelectedVariableNames = opts.SelectedVariableNames([9, 29, 31, 32, 34, 38, 39, 44, 79]);  

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
                   if (i == 1)
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
                   bomref=find(AllBomData.Date==ThisDate);
                   if isempty(bomref)
                       fprintf("Error missing day %s in BOM Data\n",ThisDate)
                       continue;
                   end
                   sampletracks=sampletracks+1;
                   minalt=min(minalt,min(T1.GeoHeight));
                   maxalt=max(maxalt,max(T1.GeoHeight));

                   col_approach=zeros(samples,1);
                   col_inb=find(strcmp(T1.App,'INB')==0);
                   col_approach(col_inb,1)=1; % 1==OUT Bound, 0==IN bound
                   %col1=zeros(samples,1);
                   %col2=ones(samples,1);
                   colour_approach=[zeros(samples,1) ones(samples,1) col_approach];

                   figure(1);
                   plot(T1.RFS_Range,T1.RangeError,'LineWidth', 2);
                   %scatter(T.RFS_Range,T.RangeError,1,colr, 'LineWidth', 2);
                   hold on;


                   figure(2);
                   coefficients = polyfit(T1.RFS_Range,T1.RangeError-RangeBias, 2);
                   xFit = linspace(min_range, max_range, samples);
                   yFit = polyval(coefficients , xFit);
                   txt=sprintf("%d",samples);
                   scatter(xFit, yFit,1,colour_approach, 'LineWidth', 2);
                   hold on;

                   figure(3);
                   xtrap = [0 50 100 150 200 250];
                   [colA altA]=getAlt(5, T1);
                   [colB altB]=getAlt(50, T1);
                   [colC altC]=getAlt(100, T1);
                   [colD altD]=getAlt(150, T1);
                   [colE altE]=getAlt(200, T1);
                   [colF altF]=getAlt(250, T1);
                   colbar=[colbar; colA; colB; colC; colD; colE; colF];
                   altbar=[altbar; altA; altB; altC; altD; altE; altF];

                   ytrap = polyval(coefficients, xtrap);
                   extrap_rng_err_0=[extrap_rng_err_0; ytrap(1)];
                   extrap_rng_err_50=[extrap_rng_err_50; ytrap(2)];
                   extrap_rng_err_100=[extrap_rng_err_100; ytrap(3)];
                   extrap_rng_err_150=[extrap_rng_err_150; ytrap(4)];
                   extrap_rng_err_200=[extrap_rng_err_200; ytrap(5)];
                   extrap_rng_err_250=[extrap_rng_err_250; ytrap(6)];
                   scatter(xtrap, ytrap, 20, [colA; colB; colC; colD; colE; colF],"filled");
                   hold on;
                   today=posixtime(datetime(ThisDate))
                   todaylessday=today-86400
                   yester=datetime(todaylessday, 'convertfrom', 'posixtime', 'Format', 'yyyy-MM-dd');

                   yesbomref=find(AllBomData.Date==yester);
                   if isempty(yesbomref)
                       fprintf("Error missing yesterday %s in BOM Data\n",yester)
                       continue;
                   end

                   if (q == 1)
                        Period="3pm "+string(yester)+" "+string(AllBomData{yesbomref,11})+"C"+string(AllBomData{yesbomref,12})+"% ";
                   else
                       Period = "9am "+ThisDate+" "+string(AllBomData{bomref,5})+"C "+string(AllBomData{bomref,6})+"% ";
                   end
                   ttrk=sprintf("%s (+/- 1 Hr)\nSample %d Aircraft\n",Period,sampletracks);
                   figure(1)
                   grid on;
                   title("Range VS Range Error : "+ttrk)
                   xlabel("Range (Nautical Miles)");
                   ylabel("Range Error (metres)")
                   %inoutcolr=[0 1 0; 0 1 1];
                   %colormap(inoutcolr);
                   %colorbar

                   figure(2)
                   title("Range VS Interpolated Range Error (Incoming/Outgoing)")
                   xlabel("Range (Nautical Miles)");
                   ylabel("Range Error (metres) (Interpolated)")
                   inoutcolr=[0 1 0; 0 1 1];
                   colormap(inoutcolr);
                   colorbar
                   grid on

                   figure(3);
                   title("Population of Range VS Interpolated Range Error at Altitude")
                   xlabel("Range (Nautical Miles)");
                   ylabel("Range Error (metres) (Interpolated)")
                   %altcolor=[[1[minalt:1:maxalt][minalt:1:maxalt]];
                   [saltbar i]=sort(altbar)
                   newcolbar=colbar(i,:)
                   colormap(newcolbar);
                   colorbar
                   set(gca,'clim',saltbar([1,end]))

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

                   figure(4)
                   x=[0 50 100 150 200 250];
                   y=[rerr_0_avg rerr_50_avg rerr_100_avg rerr_150_avg rerr_200_avg rerr_250_avg];
                   ym=[rerr_0_median rerr_50_median rerr_100_median rerr_150_median rerr_200_median rerr_250_median];
                   errorbar(x,y,[rerr_0_std rerr_50_std rerr_100_std rerr_150_std rerr_200_std rerr_250_std ]);
                   hold on
                   plot(x,ym, '-');
                   legend("Average of population at Range X","Median of population at Range X");
               end
           end
       end
    end
  
end

function [avg stdev median1]=stats(array)
   avg=mean(array);
   stdev=std(array)
   median1=median(array);
end

function  colvect=colorAtAlt(i)
   newalt = i*65535/15000
   c1=newalt/256
   c2 = rem(newalt,256)
   colvect=[1 c2/256 c1/256]
end

function [col avgalt]=getAlt(i, T)
   [TR501 fiftyidx]=find(T.RFS_Range>i-0.1);
   if isempty(fiftyidx)
      col=[1 1 1];
      avgalt=-1
   else
      [TR502 fiftyidx2]=find(TR501 < i+0.1);
      if isempty(fiftyidx2)
         col = [1 1 1];
         avgalt=-1
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

