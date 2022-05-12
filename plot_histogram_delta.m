clc
clear all
RRR={'piControl','1pctCO2'}
run1=char(RRR(1));
run2=char(RRR(2));

LAT=ncread(['/project/cmip5/hhsu/IPSL-CM6A-LR/mrsos_IPSL-CM6A-LR_' run1 '_r1i1p1f1_regrided2x2_ng.nc'],'lat');
hexMap = {'969696','fee6ce','fdd0a2','fdae6b','fd8d3c','f16913','d94801','a63603','7f2704'}
myColorMap = zeros(length(hexMap), 3); % Preallocate
for k = 1 : length(hexMap)
	thisCell = hexMap{k};
        r = hex2dec(thisCell(1:2));
        g = hex2dec(thisCell(3:4));
        b = hex2dec(thisCell(5:6));
        myColorMap(k, :) = [r, g, b];
end
ccc = myColorMap / 255; % Normalize to range 0-1



	j=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_con_full');
	j(:,1:15)=nan;
	j(:,75:90)=nan;
	AA=permute(repmat(LAT,[1 180]),[2 1]);
        for i=1:7
        b=nan(180,90);
        a=find(j(:)==i-1);
        b(a)=j(a);
        b(find(b(:)>-1))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
	for i=8
	b=nan(180,90);
        a=find(j(:)>6);
        b(a)=j(a);
        b(find(b(:)>-1))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
	allu= double(allu)
	allu=allu/sum(allu(:));
	myColor = rand(10,3); % 10 bins/colors with random r,g,b for each
	b = barh(allu, 'facecolor', 'flat');
	b.CData = ccc(1:8,:);
	b.LineWidth =2.5
	axis([0 1 0.5 9])
	set(gca, 'YTick', [1,2,8],'YTickLabel',{'0','1','\geq7'});
	set(gca, 'XTick', [0.2,0.4],'XTickLabel',{'0.2','0.4'});
	set(gca, 'box', 'off');
	set(gcf,'Units','centimeters','position',[1 1 7.2 2.8]);
        saveas(gcf,['./histogram_full.png'])



	figure
        j=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_con_wet');
        j(:,1:15)=nan;
        j(:,75:90)=nan;
        AA=permute(repmat(LAT,[1 180]),[2 1]);
        for i=1:7
        b=nan(180,90);
        a=find(j(:)==i-1);
        b(a)=j(a);
        b(find(b(:)>-1))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        for i=8
        b=nan(180,90);
        a=find(j(:)>6);
        b(a)=j(a);
        b(find(b(:)>-1))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        allu= double(allu)
        allu=allu/sum(allu(:));
        myColor = rand(10,3); % 10 bins/colors with random r,g,b for each
        b = barh(allu, 'facecolor', 'flat');
        b.CData = ccc(1:8,:);
        b.LineWidth =2.5
	axis([0 1 0.5 9])
        set(gca, 'YTick', [1,2,8],'YTickLabel',{'0','1','\geq7'});
        set(gca, 'XTick', [0.2,0.4],'XTickLabel',{'0.2','0.4'});
        set(gca, 'box', 'off');
        set(gcf,'Units','centimeters','position',[1 1 7.2 2.8]);
        saveas(gcf,['./histogram_wet.png'])

	figure
        j=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_con_tran');
        j(:,1:15)=nan;
        j(:,75:90)=nan;
        AA=permute(repmat(LAT,[1 180]),[2 1]);
        for i=1:7
        b=nan(180,90);
        a=find(j(:)==i-1);
        b(a)=j(a);
        b(find(b(:)>-1))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        for i=8
        b=nan(180,90);
        a=find(j(:)>6);
        b(a)=j(a);
        b(find(b(:)>-1))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        allu= double(allu)
        allu=allu/sum(allu(:));
        b = barh(allu, 'facecolor', 'flat');
        b.CData = ccc(1:8,:);
        b.LineWidth =2.5
        axis([0 1 0.5 9])
        set(gca, 'YTick', [1,2,8],'YTickLabel',{'0','1','\geq7'});
        set(gca, 'XTick', [0.2,0.4],'XTickLabel',{'0.2','0.4'});
        set(gca, 'box', 'off');
        set(gcf,'Units','centimeters','position',[1 1 7.2 2.8]);
        saveas(gcf,['./histogram_tran.png'])

	figure
        j=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_con_dry');
        j(:,1:15)=nan;
        j(:,75:90)=nan;
        AA=permute(repmat(LAT,[1 180]),[2 1]);
        for i=1:7
        b=nan(180,90);
        a=find(j(:)==i-1);
        b(a)=j(a);
        b(find(b(:)>-1))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        for i=8
        b=nan(180,90);
        a=find(j(:)>6);
        b(a)=j(a);
        b(find(b(:)>-1))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        allu= double(allu)
        allu=allu/sum(allu(:));
        myColor = rand(10,3); % 10 bins/colors with random r,g,b for each
        b = barh(allu, 'facecolor', 'flat');
        b.CData = ccc(1:8,:);
        b.LineWidth =2.5
        axis([0 1 0.5 9])
        set(gca, 'YTick', [1,2,8],'YTickLabel',{'0','1','\geq7'});
        set(gca, 'XTick', [0.2,0.4],'XTickLabel',{'0.2','0.4'});
        set(gca, 'box', 'off');
        set(gcf,'Units','centimeters','position',[1 1 7.2 2.8]);
        saveas(gcf,['./histogram_dry.png'])



hexMap = {'014636','276419','4d9221','7fbc41','b8e186','ffffff','fde0ef','f1b6da','de77ae','c51b7d','8e0152'}
myColorMap = zeros(length(hexMap), 3); % Preallocate
for k = 1 : length(hexMap)
        thisCell = hexMap{k};
        r = hex2dec(thisCell(1:2));
        g = hex2dec(thisCell(3:4));
        b = hex2dec(thisCell(5:6));
        myColorMap(k, :) = [r, g, b];
end
ccc = myColorMap / 255; % Normalize to range 0-1


	figure
        j=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_con_full');
	j2=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_exp_full');
	j=j2-j;
        j(:,1:15)=nan;
        j(:,75:90)=nan;
        AA=permute(repmat(LAT,[1 180]),[2 1]);
	for i=1
	b=nan(180,90);
        a=find(j(:)<-4);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
	
        for i=2:10
        b=nan(180,90);
        a=find(j(:)==i-6);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        for i=11
        b=nan(180,90);
        a=find(j(:)>4);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        allu= double(allu)
        allu=allu/sum(allu(:));
        b = barh(allu, 'facecolor', 'flat');
        b.CData = ccc(1:11,:);
        b.LineWidth =2.5
        axis([0 1 0.5 11.5])
	set(gca, 'YTick', [1,6,11],'YTickLabel',{'\leq-5','0','\geq5'});
        set(gca, 'XTick', [0.2,0.4],'XTickLabel',{'0.2','0.4'});
        set(gca, 'box', 'off');
        set(gcf,'Units','centimeters','position',[1 1 7.2 2.8]);
	saveas(gcf,['./histogram_dif_full.png'])


        figure
        j=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_con_dry');
        j2=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_exp_dry');
        j=j2-j;
        j(:,1:15)=nan;
        j(:,75:90)=nan;
        AA=permute(repmat(LAT,[1 180]),[2 1]);
        for i=1
        b=nan(180,90);
        a=find(j(:)<-4);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end

        for i=2:10
        b=nan(180,90);
        a=find(j(:)==i-6);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        for i=11
        b=nan(180,90);
        a=find(j(:)>4);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        allu= double(allu)
        allu=allu/sum(allu(:));
        b = barh(allu, 'facecolor', 'flat');
        b.CData = ccc(1:11,:);
        b.LineWidth =2.5
        axis([0 1 0.5 11.5])
        set(gca, 'YTick', [1,6,11],'YTickLabel',{'\leq-5','0','\geq5'});
        set(gca, 'XTick', [0.2,0.4],'XTickLabel',{'0.2','0.4'});
        set(gca, 'box', 'off');
        set(gcf,'Units','centimeters','position',[1 1 7.2 2.8]);
        saveas(gcf,['./histogram_dif_dry.png'])


        figure
        j=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_con_tran');
        j2=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_exp_tran');
        j=j2-j;
        j(:,1:15)=nan;
        j(:,75:90)=nan;
        AA=permute(repmat(LAT,[1 180]),[2 1]);
        for i=1
        b=nan(180,90);
        a=find(j(:)<-4);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end

        for i=2:10
        b=nan(180,90);
        a=find(j(:)==i-6);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        for i=11
        b=nan(180,90);
        a=find(j(:)>4);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        allu= double(allu)
        allu=allu/sum(allu(:));
        b = barh(allu, 'facecolor', 'flat');
        b.CData = ccc(1:11,:);
        b.LineWidth =2.5
        axis([0 1 0.5 11.5])
        set(gca, 'YTick', [1,6,11],'YTickLabel',{'\leq-5','0','\geq5'});
        set(gca, 'XTick', [0.2,0.4],'XTickLabel',{'0.2','0.4'});
        set(gca, 'box', 'off');
        set(gcf,'Units','centimeters','position',[1 1 7.2 2.8]);
        saveas(gcf,['./histogram_dif_tran.png'])


        figure
        j=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_con_wet');
        j2=ncread('/homes/hhsu/02.InfoTheo/CMIP6_RegimeShift/Analysis/Mode_Candidate.nc','Diverge_exp_wet');
        j=j2-j;
        j(:,1:15)=nan;
        j(:,75:90)=nan;
        AA=permute(repmat(LAT,[1 180]),[2 1]);
        for i=1
        b=nan(180,90);
        a=find(j(:)<-4);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end

        for i=2:10
        b=nan(180,90);
        a=find(j(:)==i-6);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        for i=11
        b=nan(180,90);
        a=find(j(:)>4);
        b(a)=j(a);
        b(find(b(:)>-99))=1;
        b=b.*(abs(cos(AA*pi/2/90)));
        allu(i)=nansum(b(:));
        end
        allu= double(allu)
        allu=allu/sum(allu(:));
        b = barh(allu, 'facecolor', 'flat');
        b.CData = ccc(1:11,:);
        b.LineWidth =2.5
        axis([0 1 0.5 11.5])
        set(gca, 'YTick', [1,6,11],'YTickLabel',{'\leq-5','0','\geq5'});
        set(gca, 'XTick', [0.2,0.4],'XTickLabel',{'0.2','0.4'});
	set(gca, 'box', 'off');
        set(gcf,'Units','centimeters','position',[1 1 7.2 2.8]);
        saveas(gcf,['./histogram_dif_wet.png'])

