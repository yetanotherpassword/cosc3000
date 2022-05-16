clear all
close all

%prefix="041_HNS";
prefix="049_MCD";

if prefix=="041_HNS";
    site="Hahn Tableland (Cairns)";
    combined="D:/Downloads/cosc3000-main/Combined/HNS_BOM_RAD_Combined/HNS_Combined.csv";
else
    site="Mt Macedon (Melbourne)";
    combined="D:/Downloads/cosc3000-main/Combined/MCD_BOM_RAD_Combined/MCD_Combined.csv";
    sic="049";
end

varNames={'Date','RangeStdDev','RangeGain','RangeBias','MinTemp','MaxTemp','am9Temp','am9RelHumid','pm3Temp','pm3RelHumid' };
opts = detectImportOptions(combined,'Delimiter',',','PartialFieldRule','fill','VariableNamingRule','preserve');
opts.SelectedVariableNames = opts.SelectedVariableNames([2, 4, 5, 6, 45, 46, 53, 54, 59, 60]);
BOMCombo = readtable(combined,opts);
BOMCombo.Properties.VariableNames=varNames;
[mint idxtm]=min(BOMCombo.MinTemp);
[maxt idxtx]=min(BOMCombo.MaxTemp);
[mint9a idxtm9a]=min(BOMCombo.am9Temp);
[maxt9a idxtx9a]=max(BOMCombo.am9Temp);
[minh9a idxhm9a]=min(BOMCombo.am9RelHumid);
[maxh9a idxhx9a]=max(BOMCombo.am9RelHumid);
[mint3p idxtm3p]=min(BOMCombo.pm3Temp);
[maxt3p idxtx3p]=max(BOMCombo.pm3Temp);
[minh3p idxhm3p]=min(BOMCombo.pm3RelHumid);
[maxh3p idxhx3p]=max(BOMCombo.pm3RelHumid);
fprintf("Yearly    Min Temp was %f C on %s\n",mint, BOMCombo.Date(idxtm))
fprintf("Yearly    Max Temp was %f C on %s\n",maxt, BOMCombo.Date(idxtx))
fprintf("Daily 9am Min Temp was %f C on %s (with RelH %f %%)\n",mint9a, BOMCombo.Date(idxtm9a),BOMCombo.am9RelHumid(idxtm9a))
fprintf("Daily 9am Max Temp was %f C on %s (with RelH %f %%)\n",maxt9a, BOMCombo.Date(idxtx9a),BOMCombo.am9RelHumid(idxtx9a))
fprintf("Daily 9am Min RelH was %f %% on %s (with Temp %f C)\n",minh9a, BOMCombo.Date(idxhm9a),BOMCombo.am9Temp(idxhm9a))
fprintf("Daily 9am Max RelH was %f %% on %s (with Temp %f C)\n",maxh9a, BOMCombo.Date(idxhx9a),BOMCombo.am9Temp(idxhx9a))

fprintf("Daily 3pm Min Temp was %f C on %s (with RelH %f %%)\n",mint3p, BOMCombo.Date(idxtm3p),BOMCombo.pm3RelHumid(idxtm3p))
fprintf("Daily 3pm Max Temp was %f C on %s (with RelH %f %%)\n",maxt3p, BOMCombo.Date(idxtx3p),BOMCombo.pm3RelHumid(idxtx3p))
fprintf("Daily 3pm Min RelH was %f %% on %s (with Temp %f C)\n",minh3p, BOMCombo.Date(idxhm3p),BOMCombo.am9Temp(idxhm3p))
fprintf("Daily 3pm Max RelH was %f %% on %s (with Temp %f C)\n",maxh3p, BOMCombo.Date(idxhx3p),BOMCombo.am9Temp(idxhx3p))



RangeStdDev=BOMCombo.RangeStdDev;
RangeGain=BOMCombo.RangeGain;
RangeBias=BOMCombo.RangeBias;
MinTemp=BOMCombo.MinTemp;
MaxTemp=BOMCombo.MaxTemp;
mat=[ RangeGain   BOMCombo.pm3Temp BOMCombo.pm3RelHumid];
[coeff,score,latent,tsquared,explained] = pca(mat);
DayDirs = dir ("D:/Downloads/cosc3000-main/SummaryData/202*");
accum_trk=[];
accum_gain=[];
maxsize=0;
for j=1:1:size(DayDirs,1)

    if DayDirs(j).name=="2021-03-17" && prefix=="041_HNS"
        % The day 2021-03-17 is missing so add null 
        accum_trk=[accum_trk; zeros(1,size(accum_trk,2))];
        accum_gain=[accum_gain; zeros(1,size(accum_gain,2))];
    else


        fname = "D:/Downloads/cosc3000-main/SummaryData/"+DayDirs(j).name+"/SUMMARY/SIC_"+prefix+"_MODE_S_summary.csv";
        %fprintf(fname)

        varNames={'Track','Date','RangeStdDev','RangeGain','RangeBias','MinRange','MaxRange','MinElev','MaxElev','Entries'};
        opts = detectImportOptions(fname,'Delimiter',',','PartialFieldRule','fill','VariableNamingRule','preserve');
        opts.SelectedVariableNames = opts.SelectedVariableNames([1, 2, 4, 5, 6, 14, 15, 16, 17, 18]);
        HNSSummary = readtable(fname,opts);
        HNSSummary.Properties.VariableNames=varNames;
        HNSSummary=HNSSummary(1:size(HNSSummary,1)-8,:);
        [ix rx]=find(HNSSummary.RangeStdDev > 1000);
        if ~isempty(ix)
            HNSSummary(ix,:)=[]
        end

        [idx0 v]=find(HNSSummary.MaxRange - HNSSummary.MinRange>200);
        HNSSummary=HNSSummary(idx0,:);

        %[idx1 v]=find(HNSSummary.MinElev >= 1);
        %HNSSummary=HNSSummary(idx1,:);

        [idx2 v]=find(HNSSummary.MaxElev <=10);
        HNSSummary=HNSSummary(idx2,:);

        [idx3 v]=find(HNSSummary.Entries >100);
        HNSSummary=HNSSummary(idx3,:);

        HNSSummary.Track = hex2dec(extractBetween(HNSSummary.Track,21,26));
        thissize = size(HNSSummary.Track,1);

        delta = thissize - maxsize ;
        if isempty(accum_trk)
            accum_trk = HNSSummary.Track';
            accum_gain = HNSSummary.RangeGain';
            maxsize = thissize;
        elseif delta < 0
            accum_trk=[accum_trk; HNSSummary.Track' zeros(-delta,1)'];
            accum_gain=[accum_gain; HNSSummary.RangeGain' zeros(-delta,1)'];
        elseif delta > 0
            tmp = [ accum_trk' ; zeros(size(accum_trk,1),delta)' ];
            accum_trk = tmp';
            accum_trk = [accum_trk; HNSSummary.Track'];
            tmp = [ accum_gain' ; zeros(size(accum_gain,1),delta)' ];
            accum_gain = tmp';
            accum_gain = [accum_gain; HNSSummary.RangeGain'];
            maxsize = thissize;
        else
            accum_trk=[accum_trk; HNSSummary.Track'];
            accum_gain=[accum_gain; HNSSummary.RangeGain'];

        end

    end
    
end

[m1 i1] = min(accum_gain);
[m2 i2] = min(m1);

rowmin = i1(i2);
colmin = i2;
minval = m2;

[m1 i1] = max(accum_gain);
[m2 i2] = max(m1);

rowmax = i1(i2);
colmax = i2;
maxval = m2;

fprintf("Largest  Range Error Gain was %f m/NM on %s by 0x%6X for %s\n",maxval, BOMCombo.Date(rowmax),accum_trk(rowmax,colmax),site)
fprintf("Smallest Range Error Gain was %f m/NM on %s by 0x%06X for %s\n",minval, BOMCombo.Date(rowmin),accum_trk(rowmin,colmin),site)

mesh(BOMCombo.Date, [1:1:size(accum_gain,2)], accum_gain')
 
accum_gain(accum_gain==0)=nan;

% mesh(BOMCombo.Date, [1:1:size(accum_gain,2)], accum_gain')



relhum=[];
for i=1:1:size(accum_gain,1)
      relhum = [relhum ones(size(accum_gain,2),1)*(BOMCombo.pm3RelHumid(i)+BOMCombo.am9RelHumid(i))/2];
end

figure
ax = surf(BOMCombo.Date', 1:1:size(accum_gain,2), accum_gain',relhum,'FaceColor','textureMap','EdgeColor','none');%,'CDataMapping','scaled')



colormap(jet(4));

hcb=colorbar;
caxis([mean(relhum,'all','omitnan')-2.5*std(relhum,[],'all','omitnan') mean(relhum,'all','omitnan')+2.5*std(relhum,[],'all','omitnan')]);
ax.LineStyle = 'none';ax.EdgeColor = 'flat';ax.FaceColor = 'flat';

titlename=sprintf("Range Gain Annual Variability for %s",site);
xlabel("(Day of Year) 01-03-2021 to 28-02-2022");
ylabel("Individual Opportunity Aircraft");
zlabel("Range Error Gain m/NM");
title(titlename);
colorTitleHandle = get(hcb,'Title');
titleString = sprintf('Relative Humidity');
set(colorTitleHandle ,'String', titleString);

