clear ; close all; clc
s=[1 2 3 7 8 9 10 11 12 13 14 15 17 18 21 23 24 25 26 27 28]; %��������
n=[4 5 6 16 19 20 22];    %��ֵ����
name={'surgery','Age','Hospital Number ','rectal temperature','pulse ','respiratory rate','temperature of extremities','peripheral pulse ','mucous membranes','capillary refill time','pain','peristalsis','abdominal distension','nasogastric tube','nasogastric reflux','nasogastric reflux PH','rectal examination - feces','abdomen','packed cell volume','total protein','abdominocentesis appearance','abdomcentesis total protein','outcome','surgical lesion','type of lesion1','type of lesion2','type of lesion3','cp_data'};
[data]=xlsread('data'); %��������

%�Ա�����ԣ�����ÿ������ȡֵ��Ƶ��
Data=cell(1,28);
for k=s      
    c=tabulate(data(:,k));  
    Data{k}=c;
end
A=Data{1,2}(:,2)~=0;
Data{1,2}=Data{1,2}(A,:); 
A=Data{1,3}(:,2)~=0;
Data{1,3}=Data{1,3}(A,:);

for k=s
    fprintf('%s:\n',name{k});
    disp(Data{1,k}(:,1:2));
end

%�����ֵ
max=max(data);     
%����Сֵ
min=min(data);      
%���ֵ
mean=nanmean(data);     
%����λ��
median=nanmedian(data);    
%�����ķ�λ��
Q1=prctile(data,25);    
%�����ķ�λ��
Q3=prctile(data,75);     
%ȱʧֵ
empty=sum(isnan(data)); 
%�����С����ֵ����λ�����ķ�λ����ȱʧֵ�ĸ���
num=[max;min;mean;median;Q1;Q3;empty];  
%num(:,s)=[];

for k=n
    fprintf('%s:\n',name{k});
    fprintf('���ֵ:%4.4f\n��Сֵ:%4.4f\n��ֵ:%4.4f\n��λ��:%4.4f\n���ķ�λ��:%4.4f\n���ķ�λ��:%4.4f\nȱʧֵ:%4.4f\n',num(1,k),num(2,k),num(3,k),num(4,k),num(5,k),num(6,k),num(7,k));  
end


for k=n
    figure(1);
    subplot(1,2,1),hist(data(:,k));
    title([name{k},'��ֱ��ͼ']);
    subplot(1,2,2),qqplot(data(:,k));
    title([name{k},'��qqͼ']);
    hold on;
    pause;
    hold off;
end

    figure(1);
    subplot(2,4,1),qqplot(data(:,4));
    title([name{4},'��qqͼ']);
    subplot(2,4,2),qqplot(data(:,5));
    title([name{5},'��qqͼ']);
    subplot(2,4,3),qqplot(data(:,6));
    title([name{6},'��qqͼ']);   
    subplot(2,4,4),qqplot(data(:,16));
    title([name{16},'��qqͼ']);
    subplot(2,4,5),qqplot(data(:,19));
    title([name{19},'��qqͼ']);
    subplot(2,4,6),qqplot(data(:,20));
    title([name{20},'��qqͼ']);
    subplot(2,4,7),qqplot(data(:,22));
    title([name{22},'��qqͼ']);
    hold on;
    pause;
    hold off;

close(figure(gcf));
%���ƺ�ͼ
% for k=n
%     figure(2);
%     boxplot(data(:,k));
%     title([name{k} ,'����ͼ']);
%     hold on;
%     pause;
%     hold off;
% end
    figure(2);
    subplot(2,4,1),boxplot(data(:,4));
    title([name{4},'����ͼ']);
    subplot(2,4,2),boxplot(data(:,5));
    title([name{5},'����ͼ']);
    subplot(2,4,3),boxplot(data(:,6));
    title([name{6},'����ͼ']);   
    subplot(2,4,4),boxplot(data(:,16));
    title([name{16},'����ͼ']);
    subplot(2,4,5),boxplot(data(:,19));
    title([name{19},'����ͼ']);
    subplot(2,4,6),boxplot(data(:,20));
    title([name{20},'����ͼ']);
    subplot(2,4,7),boxplot(data(:,22));
    title([name{22},'����ͼ']);
    hold on;
    pause;
    hold off;
close(figure(gcf));

%�޳�ȱʧֵ 
nan=data;
nan(any(isnan(nan)'),:)=[];
for k=n
    figure(3);
    subplot(2,2,1),hist(data(:,k));
    title([name{k} ,'���޳�֮ǰ��ֱ��ͼ']);
    subplot(2,2,2),hist(nan(:,k));
    title([name{k} ,'���޳�֮���ֱ��ͼ']);
    subplot(2,2,3),qqplot(data(:,k));
    title([name{k} ,'���޳�֮ǰ��qqͼ']);
    subplot(2,2,4),qqplot(nan(:,k));
    title([name{k} ,'���޳�֮���qqͼ']);     
    hold on;
    pause;
    hold off;
end
close(figure(gcf));
xlswrite('nonan.xls',nan,'sheet1','A1');
%���Ƶֵ 
hf=data;
for k=n
    h=mode(hf(:,k));
    hf1=hf(:,k);
    hf1(isnan(hf1)) = h;
    figure(4);
    subplot(2,2,1),hist(data(:,k));
    title([name{k} ,'���޸�֮ǰ��ֱ��ͼ']);
    subplot(2,2,2),hist(hf1);
    title([name{k} ,'���޸�֮���ֱ��ͼ']);
    subplot(2,2,3),qqplot(data(:,k));
    title([name{k} ,'���޸�֮ǰ��qqͼ']);
    subplot(2,2,4),qqplot(hf1);
    title([name{k} ,'���޸�֮���qqͼ']);     
    hold on;
    pause;
    hold off;
end
close(figure(gcf));
%�����
analytic_mat=data;
[m, n1] = size(analytic_mat); 
ATTRIBUTE_L = 1;
ATTRIBUTE_H = 28;
standard_line = analytic_mat(11, ATTRIBUTE_L: ATTRIBUTE_H);
% �������Ծ��󣬺�������
%CORRELATION_MAT_ATTRIBUTE ��������֮�������Ծ���
% �����к��ж������ԣ�value(i,j)������i������j������ԡ����ǶԳƾ���Ӵ��
COR_SIZE = ATTRIBUTE_H - ATTRIBUTE_L + 1; % ����Ծ���Ĵ�С
cor_mat = -ones(COR_SIZE, COR_SIZE); % ��ʼ������Ծ�������Ҫȡ�������ԣ���ʼΪ��Сֵ��-1��
for i = ATTRIBUTE_L: ATTRIBUTE_H - 1
    for j = i + 1: ATTRIBUTE_H
        merge = [analytic_mat(:, i), analytic_mat(:, j)]; % ���������ϵ�������в�����
        [NaN_line, ~] = find(isnan(merge) == 1);
        merge(NaN_line, :) = []; % ɾ������NaN�����Ա���ȷ������ϵ��
        cor_indx = i - ATTRIBUTE_L + 1;
        cor_indy = j - ATTRIBUTE_L + 1; % ����Ծ����±�
        cor_mat(cor_indx, cor_indy) = corr(merge(:, 1), merge(:, 2)); % merge�����м�ȥ��NaN�������ԣ������ϵ��
        cor_mat(cor_indy, cor_indx) = cor_mat(cor_indx, cor_indy); % �Գƾ���
    end
end
cor_mat(isnan(cor_mat)) = -1;
cor_size = size(cor_mat, 1); % �����С������������Ƿ���
for i = 1: m
    for j = ATTRIBUTE_L: ATTRIBUTE_H
        if(isnan(analytic_mat(i, j)))
            [~, index] = sort(cor_mat(j - ATTRIBUTE_L + 1, :));
            index_list = fliplr(index); % sort����fliplr��ת����ɽ��򣬵õ��ο����������ȶ��б�
            flag = 0; % ��ʶ�Ƿ�ȫ�ɹ�
            for k = 1: cor_size
                ref_attr = index_list(k); % ���ڲ�ȫ�ο�������
                if(~isnan(analytic_mat(i, ref_attr)))
                    analytic_mat(i, j) = standard_line(j - ATTRIBUTE_L + 1) / standard_line(ref_attr) * ...
                        analytic_mat(i, ref_attr + ATTRIBUTE_L - 1); % ��������ȫ���ⲻ����õķ�����
                    flag = 1;
                    break;
                end
            end
            if(flag == 0)
                disp(['Insert fail at row ', num2str(i), ' col ', num2str(j)]);
                return ;
            end
        end
    end
end
hf=analytic_mat;
for k=n
    figure(5);
    subplot(2,2,1),hist(data(:,k));
    title([name{k} ,'���޸�֮ǰ��ֱ��ͼ']);
    subplot(2,2,2),hist(hf(:,k));
    title([name{k} ,'���޸�֮���ֱ��ͼ']);
    subplot(2,2,3),qqplot(data(:,k));
    title([name{k} ,'���޸�֮ǰ��qqͼ']);
    subplot(2,2,4),qqplot(hf(:,k));
    title([name{k} ,'���޸�֮���qqͼ']);     
    hold on;
    pause;
    hold off;
end
close(figure(gcf));
xlswrite('cor.xls',analytic_mat,'sheet1','A1');

%������
%CORRELATION_MAT_ATTRIBUTE ��������֮�������Ծ���
%   �����к��ж������ԣ�value(i,j)������i������j������ԡ����ǶԳƾ���Ӵ��
analytic_mat=data;
COR_SIZE = ATTRIBUTE_H - ATTRIBUTE_L + 1; % ����Ծ���Ĵ�С

cor_mat = -ones(COR_SIZE, COR_SIZE); % ��ʼ������Ծ�������Ҫȡ�������ԣ���ʼΪ��Сֵ��-1��
for i = ATTRIBUTE_L: ATTRIBUTE_H - 1
    for j = i + 1: ATTRIBUTE_H
        merge = [analytic_mat(:, i), analytic_mat(:, j)]; % ���������ϵ�������в�����
        [NaN_line, ~] = find(isnan(merge) == 1);
        merge(NaN_line, :) = []; % ɾ������NaN�����Ա���ȷ������ϵ��
        cor_indx = i - ATTRIBUTE_L + 1;
        cor_indy = j - ATTRIBUTE_L + 1; % ����Ծ����±�
        cor_mat(cor_indx, cor_indy) = corr(merge(:, 1), merge(:, 2)); % merge�����м�ȥ��NaN�������ԣ������ϵ��
        cor_mat(cor_indy, cor_indx) = cor_mat(cor_indx, cor_indy); % �Գƾ���
    end
end
cor_mat(isnan(cor_mat)) = -1;
%SIMILARITY_MAT_SAMPLE �������������������ԡ���correlation_mat_attribute���ơ�
%   �˴�������ʵ��������ŷ����þ��룬���ԽСԽ���ơ������������������ƣ����������ע�͡�
SIM_SIZE = size(analytic_mat, 1); % ���ƾ����С����analytic_mat������һ��
sim_mat = ones(SIM_SIZE, SIM_SIZE) * 999; % ��ʼ��Ϊ������
for i = 1: SIM_SIZE - 1
    for j = i + 1: SIM_SIZE
        merge = [analytic_mat(i, ATTRIBUTE_L: ATTRIBUTE_H)', ...
            analytic_mat(j, ATTRIBUTE_L: ATTRIBUTE_H)']; % ����������ת�úϲ�Ϊ������x2�ľ���
        [NaN_line, ~] = find(isnan(merge) == 1);
        merge(NaN_line, :) = [];     
        sim_mat(i, j) = norm(merge(:, 1) - merge(:, 2)); % ��������ŷ����þ���
        sim_mat(j, i) = sim_mat(i, j); % �Գƾ���
    end
end
        sim_size = size(sim_mat, 1); % �����С������������Ƿ���
        for i = 1: m
            for j = ATTRIBUTE_L: ATTRIBUTE_H
                if(isnan(analytic_mat(i, j)))
                    [~, index_list] = sort(sim_mat(i, :));
                    flag = 0; % ��ʶ�Ƿ�ȫ�ɹ�
                    for k = 1: sim_size
                        ref_samp = index_list(k); % ���ڲ�ȫ�ο�������
                        if(~isnan(analytic_mat(ref_samp, j)))
                            analytic_mat(i, j) = analytic_mat(ref_samp, j); % ԭ�����ϣ���ȫ
                            flag = 1;
                            break;
                        end
                    end
                    if(flag == 0)
                        disp(['Insert fail at row ', num2str(i), ' col ', num2str(j)]);
                        return ;
                    end
                end
            end
        end
hf=analytic_mat;
for k=n
    figure(6);
    subplot(2,2,1),hist(data(:,k));
    title([name{k} ,'���޸�֮ǰ��ֱ��ͼ']);
    subplot(2,2,2),hist(hf(:,k));
    title([name{k} ,'���޸�֮���ֱ��ͼ']);
    subplot(2,2,3),qqplot(data(:,k));
    title([name{k} ,'���޸�֮ǰ��qqͼ']);
    subplot(2,2,4),qqplot(hf(:,k));
    title([name{k} ,'���޸�֮���qqͼ']);     
    hold on;
    pause;
    hold off;
end
close(figure(gcf));
xlswrite('sim.xls',analytic_mat,'sheet1','A1');





