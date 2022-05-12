clc
close all
clear all
addpath /homes/hhsu/Matlab_tool
FileName={'IPSL-CM6A-LR','CanESM5','AWI-ESM-1-1-LR','CMCC-ESM2','CNRM-CM6-1-HR','NorESM2-MM','MRI-ESM2-0','MIROC6'}
 load GSHHS_COAST220HL1
        GSHHS_COAST220HL2=GSHHS_COAST220HL1;
        GSHHS_COAST220HL1(:,1)=GSHHS_COAST220HL1(:,1)-180;
        GSHHS_COAST220HL2(:,1)=GSHHS_COAST220HL2(:,1)+180;
        MAP=cat(1,GSHHS_COAST220HL1,GSHHS_COAST220HL2);
        a=find(MAP(:,1)<-180);
        MAP(a,:)=[];
        a=find(MAP(:,1)>360);
        MAP(a,:)=[];
ccc=[102,194,165
252,141,98
141,160,203
231,138,195
166,216,84]/256;

ccc=[127,201,127
190,174,212
253,192,134
255,255,153
56,108,176]/256;
RRR={'piControl','1pctCO2'}
	run1=char(RRR(1));
	run2=char(RRR(2));

	ic_con_model=[18255 10950 10950 10950 10950 10950 10950 10950];
	ic_exp_model=[36518 36865 11315 18250 36517 10950 10950 10950];
for Model=5:8
	ic_con=ic_con_model(Model);
	ic_exp=ic_exp_model(Model);
	MODELNAME=char(FileName(Model))
	if exist(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run1 '_50.nc'], 'file')
        LON=ncread(['/project/cmip5/hhsu/' MODELNAME '/mrsos_' MODELNAME '_' run1 '_r1i1p1f1_regrided2x2_ng.nc'],'lon');
        LAT=ncread(['/project/cmip5/hhsu/' MODELNAME '/mrsos_' MODELNAME '_' run1 '_r1i1p1f1_regrided2x2_ng.nc'],'lat');
	fair_1=[0 2 4 4 6];
	fair_2=[0 6 12 12 18];
	fair_3=[0 10 20 20 30];
	
	a_con=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run1 '_50.nc'],'BIC');
	a_con(find(a_con(:)<-998000000))=nan;
	b_con=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run1 '_50.nc'],'Slope_1Seg');
	z=squeeze(a_con(:,:,2));
	z(find(b_con(:)<0))=nan;
	a_con(find(a_con(:)<-998000000))=nan;
	a_con(:,:,2)=z;
	x1_con=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run1 '_50.nc'],'BPx_2Seg_LHSflat');
        x2_con=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run1 '_50.nc'],'BPx_2Seg_RHSflat');
        x3_con=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run1 '_50.nc'],'BPx1_3Seg');
        x4_con=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run1 '_50.nc'],'BPx2_3Seg');

	dx=x4_con-x3_con;
        a5=a_con(:,:,5);
        a5(find(dx(:)<0))=nan;
        a_con(:,:,5)=a5;

	a_exp=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run2 '_50.nc'],'BIC');
        a_exp(find(a_exp(:)<-998000000))=nan;
	b_exp=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run2 '_50.nc'],'Slope_1Seg');
	z=squeeze(a_exp(:,:,2));
	z(find(b_exp(:)<0))=nan;
	a_exp(find(a_exp(:)<-998000000))=nan;
	a_exp(:,:,2)=z;
        x1_exp=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run2 '_50.nc'],'BPx_2Seg_LHSflat');
        x2_exp=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run2 '_50.nc'],'BPx_2Seg_RHSflat');
        x3_exp=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run2 '_50.nc'],'BPx1_3Seg');
        x4_exp=ncread(['/project/land/hhsu/03.CMIP6/BP_SMxLE_' MODELNAME '_' run2 '_50.nc'],'BPx2_3Seg');

	dx=x4_exp-x3_exp;
        a5=a_exp(:,:,5);
        a5(find(dx(:)<0))=nan;
        a_exp(:,:,5)=a5;


	latmax=size(a_con,2);
	lonmax=size(a_con,1);
	sm_mask=ncread(['/project/cmip5/hhsu/' MODELNAME '/mrsos_' MODELNAME '_piControl_r1i1p1f1_regrided2x2_ng.nc'],'mrsos',[1 1     1],[Inf Inf 1]);
        sm_mask(find(sm_mask(:)==0))=nan;
        sm_mask(152:169,75:90)=nan;
        sm_mask(:,84:90)=nan;

	mr3_con=nan(lonmax,latmax);
        mr3_exp=nan(lonmax,latmax);

	for la=1:latmax;
        for lo=1:lonmax
        	for j=1:5
                aaa_con=a_con(lo,la,j);
                BIC3_con(j)=aaa_con+fair_3(j);
		aaa_exp=a_exp(lo,la,j);
                BIC3_exp(j)=aaa_exp+fair_3(j);
        	end
                [ma3 md3]=(min(BIC3_con(:)));
		mr3_con(lo,la)=md3(1);
		[ma3 md3]=(min(BIC3_exp(:)));
                mr3_exp(lo,la)=md3(1);
        end
	end
	mr3_con(find(isnan(sm_mask(:))))=nan;
        mr3_exp(find(isnan(sm_mask(:))))=nan;
	mr3_con(sm_mask(:)<0)=nan;
        mr3_exp(sm_mask(:)<0)=nan;
        mr3_con(sm_mask(:)>10000)=nan;
        mr3_exp(sm_mask(:)>10000)=nan;
	mr3_con(sm_mask(:)==0)=nan;
        mr3_exp(sm_mask(:)==0)=nan;

figure
	subplot(2,1,1)
        a1=mr3_con(1:lonmax*0.5,:);
        a2=mr3_con(lonmax*0.5+1:lonmax,:);
        aa_con=cat(1,a2,a1);
        pcolor(LON,LAT,aa_con'); shading flat;
        hold on
        title([MODELNAME ' piControl'])
        colormap(gca,ccc)
	caxis([.5 5.5])
        hold on
	axis([0 360 -60 90])
        plot(MAP(:,1),MAP(:,2),'color',[0.7 0.7 0.7]);
	set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',15)
	set(gcf,'Units','centimeters','position',[1 1 16 14]);
        cbh = colorbar('Position', [0.175  0.05  0.65  0.025],'location','southoutside')
        cbh.Ticks = linspace(1, 5, 5) ; %Create 8 ticks from zero to 1
        cbh.TickLabels = {'001','010','110','011','111'} ;

	subplot(2,1,2)
        a1=mr3_exp(1:lonmax*0.5,:);
        a2=mr3_exp(lonmax*0.5+1:lonmax,:);
        aa_exp=cat(1,a2,a1);
        pcolor(LON,LAT,aa_exp'); shading flat;
        hold on
	axis([0 360 -60 90])
        title([MODELNAME ' 1pctCO2'])
        colormap(gca,ccc)
        caxis([.5 5.5])
        hold on
        plot(MAP(:,1),MAP(:,2),'color',[0.7 0.7 0.7]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',15)
	saveas(gcf,['./BestModel_' MODELNAME '_50.png'])
figure

	difmr3=aa_exp-aa_con;
	pcolor(LON,LAT,aa_con'); shading flat;
        hold on
	caxis([.5 5.5])
        axis([0 360 -60 90])
        plot(MAP(:,1),MAP(:,2),'color',[0.7 0.7 0.7]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',15)
        set(gcf,'Units','centimeters','position',[1 1 40 20]);
        colormap(gca,ccc)
	a=find(difmr3(:)==0);
	aa_exp(a)=0;
	[jj,pp]=find(aa_exp==1);
        jj=(LON(jj)); pp=(LAT(pp));
        plot(jj+1,pp+1,'.','MarkerSize',4,'Color',[127,201,127]/256)
	clear jj pp
	[jj,pp]=find(aa_exp==2);
        jj=(LON(jj)); pp=(LAT(pp));
        plot(jj+1,pp+1,'o','MarkerSize',4,'Color',[190,174,212]/256)
        clear jj pp
	[jj,pp]=find(aa_exp==3);
        jj=(LON(jj)); pp=(LAT(pp));
        plot(jj+1,pp+1,'x','MarkerSize',4,'Color',[253,192,134]/256)
        clear jj pp
	[jj,pp]=find(aa_exp==4);
        jj=(LON(jj)); pp=(LAT(pp));
        plot(jj+1,pp+1,'+','MarkerSize',4,'Color',[255,255,153]/256)
        clear jj pp
	[jj,pp]=find(aa_exp==5);
        jj=(LON(jj)); pp=(LAT(pp));
        plot(jj+1,pp+1,'p','MarkerSize',4,'Color',[56,108,176]/256)
        clear jj pp
	saveas(gcf,['./RegimeShift_' MODELNAME '_50.png'])

	sm_cli_con=nan(180,90);
        sm_cli_exp=nan(180,90);
	sm_cli_dif=nan(180,90);
	WT_con=nan(180,90);
	WT_exp=nan(180,90);
	CS_con=nan(180,90);
	CS_exp=nan(180,90);
	dry_con=nan(180,90);
	tran_con=nan(180,90);
	wet_con=nan(180,90);
	dry_exp=nan(180,90);
        tran_exp=nan(180,90);
        wet_exp=nan(180,90);
%
	for lat=15:90
		for lon=1:180
                sm_con=ncread(['/project/cmip5/hhsu/' MODELNAME '/mrsos_' MODELNAME '_piControl_r1i1p1f1_regrided2x2_ng.nc'],'mrsos',[lon lat     1],[1 1 1]);
                if sm_con(1)>0
	        sm_con=ncread(['/project/cmip5/hhsu/' MODELNAME '/mrsos_' MODELNAME '_piControl_r1i1p1f1_regrided2x2_ng.nc'],'mrsos',[lon lat ic_con],[1 1 18250]);
                sm_exp=ncread(['/project/cmip5/hhsu/' MODELNAME '/mrsos_' MODELNAME '_1pctCO2_r1i1p1f1_regrided2x2_ng.nc'],'mrsos',[lon lat ic_exp],[1 1 18250]);
		sm_cli_con(lon,lat)=nanmean(sm_con(:));
                sm_cli_exp(lon,lat)=nanmean(sm_exp(:));
		sm_cli_dif(lon,lat)=sm_cli_exp(lon,lat)-sm_cli_con(lon,lat);

		if mr3_con(lon,lat) == 3
                WT_con(lon,lat)=x1_con(lon,lat);
	        end
                if mr3_con(lon,lat) == 4
                CS_con(lon,lat)=x2_con(lon,lat);
                end
                if mr3_con(lon,lat) == 5
                WT_con(lon,lat)=x3_con(lon,lat);
		CS_con(lon,lat)=x4_con(lon,lat);
                end

		if mr3_exp(lon,lat) == 3
                WT_exp(lon,lat)=x1_exp(lon,lat);
                end
                if mr3_exp(lon,lat) == 4
                CS_exp(lon,lat)=x2_exp(lon,lat);
                end
                if mr3_exp(lon,lat) == 5
                WT_exp(lon,lat)=x3_exp(lon,lat);
                CS_exp(lon,lat)=x4_exp(lon,lat);
                end

		if WT_con(lon,lat)>0
		dry_con(lon,lat)=numel(find(sm_con(:)<WT_con(lon,lat)));
		else
		dry_con(lon,lat)=0;
		end
		if CS_con(lon,lat)>0
                wet_con(lon,lat)=numel(find(sm_con(:)>CS_con(lon,lat)));
		else
                wet_con(lon,lat)=0;
                end
		tran_con(lon,lat)=29218-wet_con(lon,lat)-dry_con(lon,lat);

		if WT_exp(lon,lat)>0
                dry_exp(lon,lat)=numel(find(sm_exp(:)<WT_exp(lon,lat)));
                else
                dry_exp(lon,lat)=0;
                end
                if CS_exp(lon,lat)>0
                wet_exp(lon,lat)=numel(find(sm_exp(:)>CS_exp(lon,lat)));
                else
                wet_exp(lon,lat)=0;
                end
                tran_exp(lon,lat)=29218-wet_exp(lon,lat)-dry_exp(lon,lat);
		end
		end
	end	
		save(['' MODELNAME '_diag_50_new.mat'],'sm_cli_con','sm_cli_exp','sm_cli_dif','WT_con','WT_exp','CS_con','CS_exp','dry_con','tran_con','wet_con','dry_exp','tran_exp','wet_exp');
%
end
end
