%% 清屏
clc;
clear;
close all; 

%% 参数设置
time_zero = 190250;   % 时间零点t0，可以自行设置也可以在line(32,34)读入时设置
spacing = 0.0352360817477097;                                               % 1 pixel = 0.0352360817477097 μm
filebox = ['E:\课题组资料\5-汇报类\组会\20240709-夏月星实验汇报\0708-A26\TEST'];
file = search_folder(filebox,'dat');
n_max = length(file);
px = floor(sqrt(n_max))+1; py = floor(sqrt(n_max))+1;

xl = 1;   xr = 59; yd =120; yu = 360; % xleft,xright,ydown,yup;
msdt = 8; msdp = 1;

xc = (xl+xr)/2; yc = (yd+yu)/2; % xcenter,ycenter;

%% 创建空文档，方便调用
temp_x = cell(n_max,1);
temp_y = cell(n_max,1);
temp_guass = cell(n_max,1);                                                 % 预分配一个元胞数组拟合高斯函数
tempcell = cell(1,1);

%% 命名cfsj结构组
for ii = 1:n_max
    cfsj(ii).name = cell2mat(file(ii));
    cfsj(ii).msd = 0;
    cfsj(ii).std = 0;
    cfsj(ii).deltaR = 0;
    %计算第ii个文件的延迟时间
    st_position = regexp(cfsj(ii).name,'\\');
    dot_position = regexp(cfsj(ii).name,'\.');
    decay_raw = str2double(cfsj(ii).name((st_position(end)+1):dot_position(end)-1)); % 读入的时间为6位数
    cfsj(ii).time = -(decay_raw-time_zero)/75;
end

%% 根据延迟时间给出按照顺序排列的文件名
[cfsjnewtime,id] = sort([cfsj.time],'ascend'); %id为按照ascend排序的time数组的id序列

for ii = 1:n_max
    time_tit = [num2str(cfsj(id(ii)).time) 'ps']; 
    %% 读入强度矩阵
    delta_OD = importdata(cfsj(id(ii)).name);
    [xlength,ylength] = size(delta_OD);
    %% 给出intensity矩阵
    intensity = zeros(xlength*ylength,3); 
    for yy = 1:ylength
      for xx = 1:xlength              
          intensity((yy-1).*xlength+xx,1) = xx.*spacing;
          intensity((yy-1).*xlength+xx,2) = yy.*spacing;
          intensity((yy-1).*xlength+xx,3) = delta_OD(xx,yy);
      end
    end
    %% 画图
    intensity_max = max(intensity(xlength*yd:xlength*yu,3));
    intensity_min = min(intensity(xlength*yd:xlength*yu,3));
    xx_bar = reshape(intensity(:,1),xlength,ylength)-xc*spacing;            % x轴的矩阵，y固定，x随竖轴递增
    yy_bar = reshape(intensity(:,2),xlength,ylength)-yc*spacing;            % y轴的矩阵，x固定，y随横轴递减
    subplot(px,py,ii);
    %% 在中心点附近寻找最值
    dc = 5;    % 在中心点5范围内寻找最值
    delta_R(ii) = delta_OD(xc,yc);
    sq_OD = delta_OD(xc-dc:xc+dc,yc-dc:yc+dc);
    c_pos = [xc-dc-1 yc-dc-1];
    maxr = -realmax; minr = realmax;
    for nn = 1:2*dc+1
        for mm = 1:2*dc+1
            if sq_OD(mm,nn)>maxr
                maxr = sq_OD(mm,nn);
                max_pos = [mm nn];
            end
            if sq_OD(mm,nn)<minr
                minr = sq_OD(mm,nn);
                min_pos = [mm nn];
            end
        end
    end
    if abs(minr)>abs(maxr)
        delta_R(ii) = minr;
        c_pos = min_pos+c_pos;
    else
        delta_R(ii) = maxr;
        c_pos = max_pos+c_pos;
    end
    pcolor(xx_bar(xl:xr,yd:yu),yy_bar(xl:xr,yd:yu),delta_OD(xl:xr,yd:yu));
    shading interp
    set(gcf,'Colormap',turbo);
    axis equal
    xlim([(xl-xc)*spacing (xr-xc)*spacing]);
    ylim([(yd-yc)*spacing (yu-yc)*spacing]);
    title(time_tit,'FontName','Arial','FontSize',6);
    set(gca,'CLim',[intensity_min intensity_max]);                          % 设置强度坐标

    %% 导入需要进行guass拟合的函数
    temp_guass{ii} = delta_OD(xl:xr,yd:yu);
    temp_x{ii} = xx_bar(xl:xr,yd:yu);
    temp_y{ii} = yy_bar(xl:xr,yd:yu);
    delay_time(ii) = cfsj(ii).time;                                                 % 读取对应时间坐标
end
print([filebox,'\2Dplot.tif'], '-dtiffn','-r600');

for ii = 1:n_max
    
    %% 导入拟合的xyz值，并将其平铺
    tempcell{1} = temp_guass{ii};
    zfit = cell2mat(tempcell);
    tempcell{1} = temp_x{ii};
    xfit = cell2mat(tempcell);
    tempcell{1} = temp_y{ii};
    yfit = cell2mat(tempcell);
    
    zfit = reshape(zfit,[],1);
    xfit = reshape(xfit,[],1);
    yfit = reshape(yfit,[],1);

    %% guass拟合(使用origin里的guass函数）
    A = -0.001;  xc = 0; yc = 0; z0 = 0; w = 0.03; msi = (w./2).^2;                                    % 初始参数，需要一定的调整
    ini_0 = [A xc yc z0 msi];
    
    guass_fun =@(p,q) p(4)+(p(1)./(2.*sqrt(p(5).*pi./2))).*exp(-((q(:,1)-p(2)).^2+(q(:,2)-p(3)).^2)./(2.*p(5)));
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');% 非线性最小二乘法使用levenberg-marquardt算法
    [para,resn,resi,eflag,diedai,quanzhong,jjout] = lsqcurvefit(guass_fun,ini_0,[xfit,yfit],zfit,[-20,-1,-1,-0.01,0.03],[20,1,1,0.01,1],options);
    % lsqcurvefit(函数，参数初始值，横坐标，纵坐标）返回拟合参数值;resn和resi分别代表           
    % 计算误差棒大小
    resi_var = resn./(numel(zfit)-length(ini_0));   % 样本方差 = 总方差/样本数量
    parastd = sqrt(diag(resi_var*inv(jjout'*jjout)));   % 协方差矩阵 = 样本方差*雅可比矩阵

    %% 数据处理，画图并获取msd
    A = para(1); xc = para(2); yc = para(3); z0 = para(4); msi = para(5);                    % 获取对应值
    cfsj(id(ii)).deltaR = delta_R(ii) - z0;
    zfit = zfit-z0; para(4) = 0;                      % 减去基线
    zzplot = guass_fun(para,[xfit,yfit]);
    zzplot = reshape(zzplot,[xr-xl+1,yu-yd+1]);
    zfit = reshape(zfit,[xr-xl+1,yu-yd+1]);
    xfit = reshape(xfit,[xr-xl+1,yu-yd+1]);
    yfit = reshape(yfit,[xr-xl+1,yu-yd+1]);

    onefit(ii,:) = zfit((xl+xr)./2,:);  %代表原始数据的中间列
    oneplot(ii,:) = zzplot((xl+xr)./2,:);   %代表拟合数据的中间列
    
    ttemp_file = cell2mat(file(ii));
    dot_position = regexp(ttemp_file,'\.');
    out_write = [0,yfit(1,:);xfit(:,1),zfit];
    file_out = [ttemp_file(1:dot_position(end)),'csv'];
    writematrix(out_write,file_out);
    cfsj(id(ii)).msd = msi;                                                    % 单位为μm^2
    cfsj(id(ii)).std = parastd(5);                %代表sigma的误差值
end

%% 画出x=0处的y坐标拟合图
figure
for ii = 1:n_max
    if (ii)
        subplot(px,py,ii);
        plot(yfit((xl+xr)./2,:),oneplot(ii,:),'-',yfit((xl+xr)./2,:),onefit(ii,:),'.');
        hold on
    end
end

%% 画出msd随时间的曲线
figure

%% 对start和stop位置进行选取
sr = 1; so = n_max;
for ii = 1:n_max
    if cfsj(id(ii)).time<msdt
        sr = id(ii);
    end
    if cfsj(id(ii)).time<msdp
        so = id(ii);
    end
end

%% 线性拟合y=k1*x+k2，画图
k_para = polyfit(cat(1,cfsj(sr:so).time),cat(1,cfsj(sr:so).msd),1);
errorbar(cat(1,cfsj.time),cat(1,cfsj.msd),cat(1,cfsj.std),'bo');
hold on
plot(cat(1,cfsj(sr:so).time),polyval(k_para,cat(1,cfsj(sr:so).time)),'r-');
% xlim([0 50]);

mobilty = k_para(1)*193400                                                  % 得到最后的迁移率
xlabel('Time delay (ps)','FontName','Arial','FontSize',12);             % 设置x轴标签内容和字体
ylabel('σ^2 (μm^2)','FontName','Arial','FontSize',12);             % 设置y轴标签内容和字体
xlim([cfsj(n_max).time-0.1,cfsj(1).time+0.1]);
ylim([min(cat(1,cfsj(sr:so).msd))-0.02,max(cat(1,cfsj(sr:so).msd))+0.02]);
title('MSD','FontName','Arial','FontSize',14);
text(10,min(cat(1,cfsj(sr:so).msd))-0.005,['mobility=',num2str(mobilty),'cm^2V^-^1s^-^1']);
% text(0,0,['w=',num2str(w*1000),'um']);
print([filebox,'\MSD.tif'], '-dtiffn','-r600');

%% 画图ΔR/R图
figure
plot(cat(1,cfsj.time),cat(1,cfsj.deltaR));
title('ΔR/R','FontName','Arial','FontSize',16);
xlim([cfsj(n_max).time-0.1,cfsj(1).time+0.1]);
xlabel('Time delay (ps)','FontName','Arial','FontSize',12);             % 设置x轴标签内容和字体
ylabel('ΔR/R','FontName','Arial','FontSize',12);             % 设置y轴标签内容和字体
% xlim([delay_time(1),delay_time(end)]);
print([filebox,'\deltaR.tif'], '-dtiffn','-r600');
% close all; 

%% 导出画图需要使用的数据
temp_outfile = cfsj(ii).name;
line_position = regexp(temp_outfile,'\\');
out_write = [cat(1,cfsj.time),cat(1,cfsj.deltaR),cat(1,cfsj.msd),cat(1,cfsj.std)];
file_out = [temp_outfile(1:line_position(end)),'outfigure.csv'];
writematrix(out_write,file_out);
