clc
clear all
close all

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

RRR={'piControl','1pctCO2'}
run1=char(RRR(1));
run2=char(RRR(2));
LON=ncread(['/project/cmip5/hhsu/AWI-ESM-1-1-LR/mrsos_AWI-ESM-1-1-LR_' run1 '_r1i1p1f1_regrided2x2_ng.nc'],'lon');
LAT=ncread(['/project/cmip5/hhsu/AWI-ESM-1-1-LR/mrsos_AWI-ESM-1-1-LR_' run1 '_r1i1p1f1_regrided2x2_ng.nc'],'lat');
ccc=[118,42,131
153,112,171
194,165,207
231,212,232
247,247,247
217,240,211
166,219,160
90,174,97
27,120,55]/255;
for Model=1:8
        MODELNAME=char(FileName(Model))
        load(['' MODELNAME '_diag_50_new.mat']);
	subplot(8,3,(Model-1)*3+1)
	a=(dry_exp-dry_con)/182.50;
	a1=a(1:90,:); a2=a(91:180,:); a=cat(1,a2,a1);
	pcolor(LON,LAT,a'); shading flat
	hold on
        caxis([-18 18])
        axis([45 350 -60 60])
        plot(MAP(:,1),MAP(:,2),'color',[0 0 0]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',5)
	colormap(gca,ccc)

	subplot(8,3,(Model-1)*3+2)
        a=(tran_exp-tran_con)/182.50;
        a1=a(1:90,:); a2=a(91:180,:); a=cat(1,a2,a1);
        pcolor(LON,LAT,a'); shading flat
        hold on
        caxis([-18 18])
        axis([45 350 -60 60])
        plot(MAP(:,1),MAP(:,2),'color',[0 0 0]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',5)
        colormap(gca,ccc)

	subplot(8,3,(Model-1)*3+3)
        a=(wet_exp-wet_con)/182.50;
        a1=a(1:90,:); a2=a(91:180,:); a=cat(1,a2,a1);
        pcolor(LON,LAT,a'); shading flat
        hold on
        caxis([-18 18])
        axis([45 350 -60 60])
        plot(MAP(:,1),MAP(:,2),'color',[0 0 0]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',5)
        colormap(gca,ccc)

end

        set(gcf,'Units','centimeters','position',[1 1 30 30]);
	ha=get(gcf,'children');
        a=24;
	set(ha(a),'position',[.01 .886 .32 .119])
        set(ha(a-1),'position',[.34 .886 .32 .119])
        set(ha(a-2),'position',[.67 .886 .32 .119])
        set(ha(a-3),'position',[.01 .761 .32 .119])
        set(ha(a-4),'position',[.34 .761 .32 .119])
        set(ha(a-5),'position',[.67 .761 .32 .119])
        set(ha(a-6),'position',[.01 .636 .32 .119])
        set(ha(a-7),'position',[.34 .636 .32 .119])
        set(ha(a-8),'position',[.67 .636 .32 .119])
        set(ha(a-9),'position',[.01 .511 .32 .119])
        set(ha(a-10),'position',[.34 .511 .32 .119])
        set(ha(a-11),'position',[.67 .511 .32 .119])
        set(ha(a-12),'position',[.01 .386 .32 .119])
        set(ha(a-13),'position',[.34 .386 .32 .119])
        set(ha(a-14),'position',[.67  .386 .32 .119])
        set(ha(a-15),'position',[.01  .261 .32 .119])
        set(ha(a-16),'position',[.34  .261 .32 .119])
        set(ha(a-17),'position',[.67  .261 .32 .119])
        set(ha(a-18),'position',[.01  .136 .32 .119])
        set(ha(a-19),'position',[.34  .136 .32 .119])
        set(ha(a-20),'position',[.67  .136 .32 .119])
        set(ha(a-21),'position',[.01  .001 .32 .119])
        set(ha(a-22),'position',[.34  .001 .32 .119])
        set(ha(a-23),'position',[.67  .001 .32 .119])
	saveas(gcf,['./SupFig2_shiftSM.png'])
	figure
	a=(dry_exp-dry_con)/182.50;
        a1=a(1:90,:); a2=a(91:180,:); a=cat(1,a2,a1);

	pcolor(LON,LAT,a'); shading flat
        hold on
        caxis([-18 18])
        axis([45 350 -60 60])
        plot(MAP(:,1),MAP(:,2),'color',[0 0 0]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',12)
        colormap(gca,ccc)
 	colorbar
        cbh = colorbar('Position', [0.175  0.05  0.65  0.025],'location','southoutside')
	cbh.Ticks = linspace(-18, 18, 10) ; %Create 8 ticks from zero to 1

	saveas(gcf,['./SupFig2_shiftS_colorbarM.png'])
