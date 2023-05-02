clear all
close all
clc

set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

R=@(x) arrayfun(@(x) [cosd(x), sind(x); -sind(x), cosd(x)],x, 'uni',false);
a{2}=[3/2,3/2;sqrt(3)/2,-sqrt(3)/2]; %assigns primitive for Hexagonal lattice consistent with graphene
%a{2}=a{1};
%a{3}=a{1};
%a{4}=a{1};
%a{5}=a{1};
%a{1}=[2, 0; 1, sqrt(3)]';
a{1}=eye(2);
%q{1}=rvec([1,2],3,120);
%q{2}=rvec([1,-1],sqrt(3),90);
%a{2}=[1,2;q{1}];
%a{3}=[1,-1;q{2}];
%r1=6;
%s={'6',r1,0,0,30;'6',3*r1,0,0,30;'6',6*r1,0,0,30;'6',9*r1,0,0,30;'6',12*r1,0,0,30}; %sets shapes for embedding lattices
w=cell([1,size(a{1},1)]);
[w{:}]=deal(-3:3);
L{1}=lat(a,w);
%L{2}=lat({reciprocal(a{1})},w);
%r={0:30,0:60,0:120,0:180};
 %f=@(x,y) cos(2*pi()*x)^2*sin(2*pi()*y)^2;
 %F=arrayfun(@(x,y) f(x,y),L{1}(1,:,1),L{1}(2,:,1));
% G=arrayfun(@(x,y) f(x,y),L{2}(1,:,1),L{2}(2,:,1));
%k=sqrt(size(F,2));
fig1=figure;
%surf(reshape(L{1}(1,:,1),k,[]),reshape(L{1}(2,:,1),k,[]),reshape(F,k,[]),'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud');
%hold on
scatter(L{1}(1,:,1),L{1}(2,:,1), 7.5,'filled','DisplayName','$\mathcal{R}_1$');
%exportgraphics(fig1,'latfun1.pdf','BackgroundColor','none');
% fig2=figure;
% surf(reshape(L{2}(1,:,1),k,[]),reshape(L{2}(2,:,1),k,[]),reshape(G,k,[]),'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud');
% hold on
% scatter(L{2}(1,:,1),L{2}(2,:,1), 7.5,'filled','DisplayName','$\mathcal{R}_1$');
% exportgraphics(fig2,'latfun2.pdf','BackgroundColor','none');

%tplot(L,r);
% L{2}=twist(L{1},r);
% fig=figure;
% scatter(L{1}(1,:,1),L{1}(2,:,1), 7.5,'filled','DisplayName','$\mathcal{R}_1$');
 axis off;
 set(gca, 'color', 'none');
 axis([-20 20 -20 20])
% hold on
% scatter(L{2}(1,:,1),L{2}(2,:,1), 7.5,'filled','DisplayName','$\mathcal{R}_1$');
% hold off
% exportgraphics(fig,'2g9deg.pdf','BackgroundColor','none');

%for i=1:size(L{1},2)
%    f(i)=sin(2*pi*L{1}(1,i,1)/norm(a{1}(:,1)))*cos(2*pi*L{1}(2,i,1)/norm(a{1}(:,2)));
%end
%fig2=figure;
%scatter(f(1,:),f(2,:))
%surf(real(reshape(L{1}(1,:,1),13,[])),real(reshape(L{1}(2,:,1),13,[])),real(reshape(f,13,[])))



function tplot(L,r)
F=figure;
p{1}=scatter(L{1}(1,:,1),L{1}(2,:,1),'.','DisplayName','$\mathcal{R}_1$ @ $0^\circ$');
hold on

axis square;
axis equal;
axis off;
lgd=legend;
lgd.ItemTokenSize(1) = 10;
axis([-25 25 -25 25])

V=VideoWriter('lrot','MPEG-4');
V.Quality=100;
open(V);
for i=2:size(r,2)
    L{i}=twist(L{1}(:,:,i),r{:,i-1});
    T{1,i-1}=L{i}(1,:,1);
    T{2,i-1}=L{i}(2,:,1);
    p{i}=scatter(T{1,i-1},T{2,i-1},'.');
    p{i}.XDataSource='T{1,i-1}';
    p{i}.YDataSource='T{2,i-1}';
    for j=1:size(r{:,i-1},2)
        T{1,i-1}=L{i}(1,:,j);
        T{2,i-1}=L{i}(2,:,j);
        p{i}.DisplayName=strcat('$\mathcal{R}_',num2str(i),'$ @ $',num2str(r{:,i-1}(j)),'^\circ$');
        refreshdata(p{i})
        drawnow
        %pause()
        %exportgraphics(F,'lrot.gif','Append',true);
        M=getframe(F);
        writeVideo(V,M);
    end
end
%exportgraphics(F,'quadlat.pdf','BackgroundColor','none')
close(V);


end

% function rnew=crot(r) %modify to add nd rotations
%     R=@(x) arrayfun(@(x) [cosd(x), sind(x); -sind(x), cosd(x)],x, 'uni',false);
%     rnew=R(r);
%     rnew=cat(3,rnew{:});
% end
% end

% function Lnew=pformat(L) %modifys multiarray to use scatter more efficiently
% for i=1:size(L,1)
%     Lnew{i}=squeeze(L(i,:,:));
% end
% end

function a=reciprocal(a) %gives the dual primitive
    a=2*pi*inv(a)';
end
 
function w=rvec(v,n,t) %given a vector, norm, and angle, calculates the corresponding vector
R=@(x) [cosd(x), sind(x); -sind(x), cosd(x)];
if n==0
    w=(R(t)*v')';
else
    w=(n/norm(v)*R(t)*v')';
end
end

function Lnew=twist(L,r) %produces twisted lattices with input {lattices, rotation, display lattice}
    R=@(x) arrayfun(@(x) [cosd(x), sind(x); -sind(x), cosd(x)],x, 'uni',false);
    rnew=R(r);
    rnew=cat(3,rnew{:});
    Lnew=pagemtimes(rnew,L);
end

function [L]=radbed(a,n,s,t1,t2) %produces embedded lattice structure with input {primitives, window, shapes, display lattice, w/ shapes or w/o}
w=cell([1,size(a{1},1)]);
[w{:}]=deal(n:-n);
L=lat(a,w,0);
p=cell([1,size(s,1)]);
[p{:}]=deal(NaN);
k=strcmp('circle',s(:,1));
if k==0
    [p(~k)]=multipoly(s(~k,:));
elseif ~k==0
    [p(k)]=circ(s(k,2:size(s(k,:),2)));
else
    [p(k)]=circ(s(k,2:size(s(k,:),2)));
    [p(~k)]=multipoly(s(~k,:));
end
pn=cellfun(@(y) cellfun(@(x) isinterior(y,x'), L, 'uni', 0),p, 'uni', 0);
%t=mod(~pn{1,1}{1,3}+~pn{1,2}{1,3},2);
L{1}=L{1}(:,pn{1,1}{1,1});
for i=2:size(a,2)
    L{i}=L{i}(:,~pn{1,i-1}{1,i} & pn{1,i}{1,i});
end
if t1==1
    for i=1:size(L,2)
        hold on
        scatter(L{i}(1,:),L{i}(2,:),'.')
        if t2==1
            plot(p{i})
        end
    end
end
L=horzcat([L{:}]);
end

function p=multipoly(o) %produces polygons with input {#sides, radius, x-center, y-center, rotation}
for i=1:size(o,1)
    p{i}=rotate(nsidedpoly(str2num(o{i,1}),'Center',[o{i,3},o{i,4}],'Radius',o{i,2}),o{i,5});
end
end

function p=circ(o) %produces circles with input {radius, x-center, y-center, sample size}
theta = cellfun(@(x) (0:x-1)*(2*pi/x), o(:,4), 'uni', 0);
x1=cellfun(@plus, o(:,2), cellfun(@times, o(:,1),cellfun(@cos, theta, 'uni', 0), 'uni', 0), 'uni',0)';
x2=cellfun(@plus, o(:,3), cellfun(@times, o(:,1),cellfun(@sin, theta, 'uni', 0), 'uni', 0), 'uni',0)';
for i=1:size(o,1)
    p{i}=polyshape(x1{1,i},x2{1,i});
end
end

function L=lat(a,n) %produces lattices with input {primitives, window}
X=cell(size(n));
[X{:}]=ndgrid(n{:});
X=cellfun(@(x) x(:)', X, 'uni', 0);
X=vertcat(X{:});
for i=1:size(a,2)
    L(:,:,i)=a{i}*X;
end
end