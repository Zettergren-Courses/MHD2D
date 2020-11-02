%GET READY FOR THIS SET OF PLOTS
clear,clc,close all;
!rm ./plots/*.png


%OPEN FILE AND GET THE GRID DATA
fid=fopen('MHD2D.dat');
lx1=fread(fid,1,'integer*4');
x1=fread(fid,lx1+4,'real*8');
lx2=fread(fid,1,'integer*4');
x2=fread(fid,lx2+4,'real*8');

indsx1=1:20:lx1;
indsx2=1:20:lx2;
[X2,X1]=meshgrid(x2(indsx2),x1(indsx1));
[X2all,X1all]=meshgrid(x2,x1);


it=1;
while ~feof(fid)
    t=fread(fid,1,'real*8');
    if feof(fid)
      break;
    end

    n=fread(fid,(lx1+4)*(lx2+4),'real*8');
    n=reshape(n,[lx1+4,lx2+4]);
    
    v1=fread(fid,(lx1+4)*(lx2+4),'real*8');
    v1=reshape(v1,[lx1+4,lx2+4]);
        
    v2=fread(fid,(lx1+4)*(lx2+4),'real*8');
    v2=reshape(v2,[lx1+4,lx2+4]);
    
    v3=fread(fid,(lx1+4)*(lx2+4),'real*8');
    v3=reshape(v3,[lx1+4,lx2+4]);
    
    B1=fread(fid,(lx1+4)*(lx2+4),'real*8');
    B1=reshape(B1,[lx1+4,lx2+4]);
    
    B2=fread(fid,(lx1+4)*(lx2+4),'real*8');
    B2=reshape(B2,[lx1+4,lx2+4]);    
    
    B3=fread(fid,(lx1+4)*(lx2+4),'real*8');
    B3=reshape(B3,[lx1+4,lx2+4]);    
    
    p=fread(fid,(lx1+4)*(lx2+4),'real*8');
    p=reshape(p,[lx1+4,lx2+4]);
    
    T=p./n/1.38e-23;

%     g1=fread(fid,(lx1)*(lx2),'real*8');
%     g1=reshape(g1,[lx1,lx2]);
    
    
    figure(1);
    clf;
    
%     pcolor(x2/1e3,x1/1e3,p);
%     shading flat;
%     colorbar
%     xlabel('x_2 [km]');
%     ylabel('x_1 [km]');
%     hold on;
%     quiver(X2/1e3,X1/1e3,B2(indsx1,indsx2)*1e9,B1(indsx1,indsx2)*1e9,'Color','white');
%     hold off;
%     title(sprintf('pressure [Pa], B1,B2 arrows, t=%f',t))
    
%     pcolor(x2(1:lx2)/1e3,x1(1:lx1)/1e3,v3(1:lx1,1:lx2));
%     shading flat;
% %    colormap(jet(256));
%     colorbar
%     xlabel('x_2 [km]');
%     ylabel('x_1 [km]');
%     hold on;
%     %quiver(X1/1e3,X2/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
%     quiver(X2/1e3,X1/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
%     hold off;
%     title(sprintf('tracer [arb. units], v1,v2 arrows, t=%f',t));

%     pcolor(x2/1e3,x1/1e3,B3*1e9);
%     shading flat;
%     colorbar;
%     xlabel('x_2 [km]');
%     ylabel('x_1 [km]');
%     hold on;
%     quiver(X2/1e3,X1/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
%     hold off;
%     title(sprintf('B_3 [nT], v1,v2 arrows, t=%f',t));
    
%     pcolor(x2/1e3,x1/1e3,p);
%     shading flat;
%     %colormap(jet(256));
%     colorbar
%     xlabel('x_2 [km]');
%     ylabel('x_1 [km]');
%     hold on;
% %    quiver(X1/1e3,X2/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
%     quiver(X2/1e3,X1/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
%     hold off;
%     title(sprintf('pressure [Pa], v1,v2 arrows, t=%f',t));

%     T=2/3*p./n/1.38e-23;
%     pcolor(x2/1e3,x1/1e3,v1);
%     shading flat;
%     %colormap(jet(256));
%     colorbar
%     xlabel('x_2 [km]');
%     ylabel('x_1 [km]');
%     hold on;
% %    quiver(X1/1e3,X2/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
%     quiver(X2/1e3,X1/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
%     hold off;
%     title(sprintf('density [1/m^3], v1,v2 arrows, t=%f',t));
    
    
    subplot(221)
%    subplot(121);
    pcolor(x2/1e3,x1/1e3,n);
    shading flat;
    colorbar
    xlabel('x_2 [km]');
    ylabel('x_1 [km]');
    hold on;
    quiver(X2/1e3,X1/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
    hold off;
    title(sprintf('density [1/m^3], v1,v2 arrows, t=%8.2f',t));

%     subplot(221)
%     pcolor(x2/1e3,x1/1e3,log10(n));
%     shading flat;
%     colorbar
%     xlabel('x_2 [km]');
%     ylabel('x_1 [km]');
%     hold on;
%     quiver(X2/1e3,X1/1e3,v2(indsx1,indsx2)/1e3,v1(indsx1,indsx2)/1e3,'Color','white');
%     hold off;
%     title(sprintf('density [1/m^3], v1,v2 arrows, t=%8.2f',t));

    subplot(222)
%    subplot(122)
    pcolor(x2/1e3,x1/1e3,p);
    shading flat;
    colorbar
    xlabel('x_2 [km]');
    ylabel('x_1 [km]');
    hold on;
    quiver(X2/1e3,X1/1e3,B2(indsx1,indsx2)*1e9,B1(indsx1,indsx2)*1e9,'Color','white');
%     startx2=linspace(min(x2/1e3),max(x2/1e3),20);
%     startx1=min(x1/1e3)*ones(1,20);
%     h=streamline(X2all/1e3,X1all/1e3,B2,B1,startx2,startx1);
%     for il=1:numel(h)
%        h(il).Color=[1 1 1]; 
%     end
    hold off;
    title('pressure [Pa], B1,B2 arrows');

%     subplot(222)
%     pcolor(x2/1e3,x1/1e3,log10(p));
%     shading flat;
%     colorbar
%     xlabel('x_2 [km]');
%     ylabel('x_1 [km]');
%     hold on;
%     quiver(X2/1e3,X1/1e3,B2(indsx1,indsx2)*1e9,B1(indsx1,indsx2)*1e9,'Color','white');
%     hold off;
%     title('pressure [Pa], B1,B2 arrows');
    
%     subplot(223)
%     pcolor(x2/1e3,x1/1e3,v3);
%     shading flat;
%     colorbar;
%     xlabel('x_2 [km]');
%     ylabel('x_1 [km]');
%     title('u_3 [km/s]');

    subplot(223)
    pcolor(x2/1e3,x1/1e3,T);
    shading flat;
    colorbar;
    xlabel('x_2 [km]');
    ylabel('x_1 [km]');
    title('T [K]');
    
    subplot(224)
    pcolor(x2/1e3,x1/1e3,B3*1e9);
    shading flat;
    colorbar;
    xlabel('x_2 [km]');
    ylabel('x_1 [km]');
    title('B_3 [nT]')

    ndigit=log10(floor(it))+1;
    padlevel=6;
    strfile=num2str(floor(it));
    while (numel(strfile)<padlevel)
      strfile=['0',strfile];
    end
    print('-dpng',['./plots/',strfile,'.png'],'-r300');
    it=it+1;
end
fclose(fid);


% %MAKE A MOVIE WITH MPLAYER (DEFAULT MOVIE.AVI FILENAME AND 10 FPS)
% fps=10;
% currloc=pwd;
% plotloc=[currloc,'/plots'];
% moviefile=[plotloc,'/movie.avi'];
% if ismac     %OSX does weird stuff with calls to bash, could never get this to work
%   path1 = getenv('PATH');
%   path1 = [path1 ':/usr/local/bin'];
%   setenv('PATH', path1);
%   setenv('DYLD_LIBRARY_PATH', '/usr/local/bin')
% end
% system(['sh ',plotloc,'/vid_encode2.sh ',plotloc,' ',moviefile,' ',fps]);
