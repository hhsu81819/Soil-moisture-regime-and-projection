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
ccc=[178,24,43
214,96,77
244,165,130
253,219,199
247,247,247
209,229,240
146,197,222
67,147,195
33,102,172]/255;

%ccc=flipud(ccc);

aall=nan(8,180,90);
ball=nan(8,180,90);
call=nan(8,180,90);
for Model=1:8
        MODELNAME=char(FileName(Model))
        load(['' MODELNAME '_diag_50_new.mat']);
	a=(sm_cli_dif./sm_cli_con)*100;
	b=(WT_exp-WT_con)./WT_con*100;
	c=(CS_exp-CS_con)./CS_con*100;
	aall(Model,:,:)=a;
	ball(Model,:,:)=b;
	call(Model,:,:)=c;
end
	
	aall=squeeze(nanmean(aall,1));
        ball=squeeze(nanmean(ball,1));
        call=squeeze(nanmean(call,1));
	
	subplot(3,1,1)
	a1=aall(1:90,:); a2=aall(91:180,:); a=cat(1,a2,a1);
	pcolor(LON,LAT,a'); shading flat
	hold on
        caxis([-18 18])
        axis([45 350 -60 60])
        plot(MAP(:,1),MAP(:,2),'color',[0 0 0]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',10)
	colormap(gca,ccc)
	title('SM climatology change (%)')
	cbh=colorbar
	cbh.Ticks = linspace(-18, 18, 10)
	subplot(3,1,2)
        a1=ball(1:90,:); a2=ball(91:180,:); a=cat(1,a2,a1);
        pcolor(LON,LAT,a'); shading flat
        hold on
        caxis([-18 18])
        axis([45 350 -60 60])
        plot(MAP(:,1),MAP(:,2),'color',[0 0 0]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',10)
        colormap(gca,ccc)
        title('WP change (%)')
	cbh=colorbar
	cbh.Ticks = linspace(-18, 18, 10)
	subplot(3,1,3)
        a1=call(1:90,:); a2=call(91:180,:); a=cat(1,a2,a1);
        pcolor(LON,LAT,a'); shading flat
        hold on
        caxis([-18 18])
        axis([45 350 -60 60])
        plot(MAP(:,1),MAP(:,2),'color',[0 0 0]);
        set(gca,'XTick',[]);  set(gca,'YTick',[]);      set(gca,'Fontsize',10)
        colormap(gca,ccc)
        title('CSM change (%)')
	cbh=colorbar
	cbh.Ticks = linspace(-18, 18, 10)
        set(gcf,'Units','centimeters','position',[1 1 20 22]);

        saveas(gcf,['./SupFig_SM_WP_CSM.png'])
	
