clear all
close all
AllDays=dir ("~/bn_data/")
%  Z=1:1:100;
% %Colormap is defined as a 3 column matrix, each row being an RGB triplet
% map = zeros(numel(Z),3);
% map(:,1)=1;
% map(:,2)=0;
% map(:,3)=Z./max(Z);
% %Set the current Colormap
% colormap(map);
% %Display Colorbar
% colorbar
inoutcolr=[0 1 0; 0 1 1]
colormap(inoutcolr);
colorbar

altcolr=[0 1 0]

a=0
for i = 1:1:size(AllDays,1)
    if (AllDays(i).isdir && length(AllDays(i).name) == 10 && AllDays(i).name=="2021-06-06")
       ThisDate=sprintf ("%s",AllDays(i).name);
       CurrPath="~/bn_data/"+ThisDate+"/SIC_045_BNE/TGT/CSV/";
       AllTracks=dir (CurrPath);
       NumberOfTracksThisDay=size(AllTracks,1)
       for j=1:1:NumberOfTracksThisDay
           if length(AllTracks(j).name)>2
               % (>2 means: ignore current and parent dir filenames)
               ThisTrackFname=CurrPath+AllTracks(j).name
               details=split(AllTracks(j).name,"_")
               if length(details)>8
                   qos="yes"
               else
                   qos="no"
               end
               msid=string(details(6))
               varNames = {'RFS_Detection_Time','RFS_Range','RFS_Azimuth','RFS_Velocity','RFS_Heading','RangeError','App'} ;
               varTypes = {'double','double','double','double','double','double','char'} ;
               delimiter = ',';
               dataStartLine = 1;
               extraColRule = 'ignore';
               opts = detectImportOptions(ThisTrackFname,'Delimiter',delimiter,'PartialFieldRule','keep')
               opts.SelectedVariableNames = opts.SelectedVariableNames([29, 31, 32, 38, 39, 44, 79]);  

               %opts.VariableNames=varNames;
               %opts.VariableTypes=varTypes;
              % opts = delimitedTextImportOptions( 'Delimiter',delimiter,...
               %                 'DataLines', dataStartLine,...
               %        
               fid = fopen(ThisTrackFname,'rt');
summary1 = fscanf(fid, '%s\n');
 % is 
fclose(fid);
a=split(summary1,',');
RangeBias=str2double(a(163)); % ie field 71 on 2nd line
         %'ExtraColumnsRule',extraColRule); 
               
               T = readtable(ThisTrackFname,opts);

               T.Properties.VariableNames=varNames;
               %figure
               %plot(T.RFS_Range,T.RangeError)
               hold on

               coefficients = polyfit(T.RFS_Range,T.RangeError-RangeBias, 2);
               samples=size(T.RFS_Range,1);
               max_range=max(T.RFS_Range);
               min_range=min(T.RFS_Range);
               if (samples < 200 || max_range-min_range < 200)
                   continue;
               end
xFit = linspace(min_range, max_range, samples);

yFit = polyval(coefficients , xFit);
txt=sprintf("%d",samples)
col3=zeros(samples,1)
cols=find(strcmp(T.App,'INB')==0)
col3(cols,1)=1;
col1=zeros(samples,1)
col2=ones(samples,1)
colr=[col1 col2 col3]
%plot(T.RFS_Range, 'b.', 'MarkerSize', 15); 
%hold on;
%h=plot(xFit, yFit, 'LineWidth', 2); 
h=scatter(xFit, yFit,1,colr, 'LineWidth', 2);

%Colormap is defined as a 3 column matrix, each row being an RGB triplet

%Set the current Colormap

grid on;
               %plot(T.RFS_Range(find(T.App=="INB")),T.RangeError(find(T.App=="INB")))
               %plot(T.RFS_Range(find(T.App=="OUTB")),T.RangeError(find(T.App=="OUTB")))
               
           end
       end

    end
  
end

